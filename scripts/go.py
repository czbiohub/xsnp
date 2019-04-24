#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
from collections import defaultdict

ARG_MIN_DEPTH = 2
ARG_MIN_GENOME_COVERED_BASES = 10
ARG_MAX_SITE_RATIO = 10

# Thread-safe and timestamped prints.
tslock = multiprocessing.RLock()


def timestamp(t):
    # We do not use "{:.3f}".format(time.time()) because its result may be
    # up to 0.5 ms in the future due to rounding.  We want flooring here.
    s = str(int(t * 10))
    return s[:-1] + "." + s[-1:]


def tsfmt(msg):
    ts = timestamp(time.time()) + " "
    msg = ts + msg.replace("\n", "\n" + ts)
    return msg


def tsout(msg):
    with tslock:
        sys.stdout.write(str(msg))
        sys.stdout.write("\n")


def tserr(msg):
    with tslock:
        sys.stderr.write(str(msg))
        sys.stderr.write("\n")


def tsprint(msg):
    tserr(tsfmt(msg))


schema = {
   "count_a": (int, "A"),
   "count_c": (int, "C"),
   "count_g": (int, "G"),
   "count_t": (int, "T"),
   "depth":   (int,),
   "ref_pos": (int,),
   "ref_id": (lambda ref_id: ref_id.replace("|", "_"),),
   "ref_allele": (str,)
}

# A simple rows -> hashes converter.
# Credit: github.com/snayfatch/MIDAS
def parse_table(rows, schema={}):
    raw_headers = next(rows)  # pylint: disable=stop-iteration-return
    headers = []
    functions = []
    for rh in raw_headers:
        f = lambda x: x
        h = rh
        if rh in schema:
            dt = schema[rh]
            if len(dt) > 0:
                f = dt[0]
            if len(dt) > 1:
                h = dt[1]
        functions.append(f)
        headers.append(h)
    yield dict(zip(headers, range(len(headers))))
    for values in rows:
        assert len(headers) == len(values)
        yield [f(v) for f,v in zip(functions, values)]

def tsv_rows(path):
    # TODO:  Support s3 and compressed files.
    with open(path, "r") as stream:
        for line in stream:
            yield line.rstrip("\n").split("\t")


def print_top(counters, how_many=5):
    print(json.dumps(sorted(((depth, contig_id) for contig_id, depth in counters.items()), reverse=True)[:how_many], indent=4))

def chomp(s, ending):
    assert s.endswith(ending)
    return s[:-len(ending)]

def process(sample_pileup_path):
    tsprint(f"Processing sample pileup path: {sample_pileup_path}")
    sample_name = chomp(sample_pileup_path, ".pileup")
    paramstr = f"gcb{ARG_MIN_GENOME_COVERED_BASES}.dp{ARG_MIN_DEPTH}.sr{ARG_MAX_SITE_RATIO}"
    output_path_site_id = f"{sample_name}.sites.{paramstr}.tsv"
    output_path_contig_stats = f"{sample_name}.contig_stats.tsv"
    t_start = time.time()
    sites = {}
    contig_depth = defaultdict(int)
    genome_covered_bases = defaultdict(int)
    contig_covered_bases  = defaultdict(int)
    table_iterator = parse_table(tsv_rows("SRS023176.pileup"), schema)
    columns = next(table_iterator)
    # Add derived columns
    columns["site_id"] = len(columns)
    columns["genome_id"] = len(columns)
    columns["number_alleles"] = len(columns)
    # Get integer keys for columns
    c_ref_id = columns["ref_id"]
    c_depth = columns["depth"]
    c_ref_pos = columns["ref_pos"]
    c_ref_allele = columns["ref_allele"]
    c_A = columns["A"]
    c_C = columns["C"]
    c_G = columns["G"]
    c_T = columns["T"]
    c_site_id = columns["site_id"]
    c_genome_id = columns["genome_id"]
    c_number_alleles = columns["number_alleles"]
    for line, row in enumerate(table_iterator):
        if line % (1000*1000) == 0:
            tsprint(f"Processing {line}.")
        # if line == (2*1000*1000):
        #    break
        contig_id = row[c_ref_id]
        row_depth = row[c_depth]
        # Add derived columns
        site_id = contig_id + "|" + str(row[c_ref_pos]) + "|" + str(row[c_ref_allele])
        genome_id = contig_id.split("_", 1)[0]
        assert len(row) == c_site_id
        row.append(site_id)
        assert len(row) == c_genome_id
        row.append(genome_id)
        number_alleles = 0
        for nt_count in (c_A, c_C, c_G, c_T):
            if row[nt_count] / row_depth >= 0.01:
                number_alleles += 1
        assert len(row) == c_number_alleles
        row.append(number_alleles)
        # Index and aggregate rows
        sites[site_id] = row
        contig_depth[contig_id] += row_depth
        contig_covered_bases[contig_id] += 1
        genome_covered_bases[genome_id] += 1
        if line < 10:
            print(json.dumps(dict(zip(columns.keys(), row)), indent=4))
    # print_top(contig_depth)
    # print_top(contig_covered_bases)
    # print_top(genome_covered_bases)
    with open(output_path_site_id, "w") as o1:
        o1.write("site_id\n")
        output_sites = 0
        for site_id, row in sites.items():
            if row[c_depth] < ARG_MIN_DEPTH:
                continue
            if genome_covered_bases[row[c_genome_id]] < ARG_MIN_GENOME_COVERED_BASES:
                continue
            contig_id = row[c_ref_id]
            average_contig_depth = contig_depth[contig_id] / contig_covered_bases[contig_id]
            site_ratio = row[c_depth] / average_contig_depth
            if site_ratio > ARG_MAX_SITE_RATIO:
                continue        
            o1.write(site_id + "\n")
            output_sites += 1
    with open(output_path_contig_stats, "w") as o2:
        o2.write("contig_id\tgenome_id\tcontig_total_depth\tcontig_covered_bases\n")
        for contig_id, contig_total_depth in contig_depth.items():
            covered_bases = contig_covered_bases[contig_id]
            genome_id = contig_id.split("_", 1)[0]
            o2.write(f"{contig_id}\t{genome_id}\t{contig_total_depth}\t{covered_bases}\n")
    t_end = time.time()
    tsprint(f"Output {output_sites} sites passing filters, out of {len(sites)} total sites.  Pass rate: {output_sites/len(sites)*100:3.1f} percent.")
    tsprint(f"Run time {t_end - t_start} seconds, or {len(sites)/(t_end - t_start):.1f} sites per second.")


if __name__ == "__main__":
    assert len(sys.argv) > 1
    for sample_pileup_path in sys.argv[1:]:
        process(sample_pileup_path)
