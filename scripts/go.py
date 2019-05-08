#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
from collections import defaultdict


import param
from util import *


def process(sample_pileup_path, num_threads, thread_id, contig_accumulator, genome_accumulator):

    #tsprint(f"{sample_pileup_path}_tid{thread_id}: Processing sample pileup path")
    sample_name = chomp(sample_pileup_path, ".pileup")
    paramstr = f"gcb{param.MIN_GENOME_COVERED_BASES}.dp{param.MIN_DEPTH}.band{thread_id}"
    banded_output_path = f"banded/{sample_name}.sites.{paramstr}.tsv"
    t_start = time.time()

    sites = {}
    sites_count = 0

    table_iterator = parse_table(tsv_rows_slice(sample_pileup_path, num_threads, thread_id), param.sample_pileup_schema)
    columns = next(table_iterator)

    # Add derived columns
    columns["site_id"] = len(columns)
    columns["genome_id"] = len(columns)
    columns["number_alleles"] = len(columns)
    columns['nz_allele'] = len(columns)
    columns['nz_allele_count'] = len(columns)

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
    c_nz_allele = columns["nz_allele"]
    c_nz_allele_count = columns["nz_allele_count"]

    # Ouput column indices
    oc_contig_depth, oc_contig_covered_bases = range(2)
    og_genome_depth, og_genome_covered_bases = range(2)

    contig_acc = contig_accumulator[sample_name]
    genome_acc = genome_accumulator[sample_name]

    for line, row in enumerate(table_iterator):
        if line % (1000*100) == 0:
            tsprint(f"{sample_pileup_path}_tid{thread_id}: Processing {line}.")
        if line == param.MAX_LINES:
            break

        # Unpack frequently accessed columns.
        contig_id = row[c_ref_id]
        depth = row[c_depth]
        ref_pos = row[c_ref_pos]
        ref_allele = row[c_ref_allele]

        # Compute derived columns.
        site_id = f"{contig_id}|{ref_pos}|{ref_allele}"
        genome_id = contig_id.split("_", 1)[0]

        # Index and aggregate rows
        sites_count += 1

        acc1 = contig_acc.get(contig_id)
        if acc1:
            acc1[oc_contig_depth] += depth
            acc1[oc_contig_covered_bases] += 1
        else:
            acc1 = [depth, 1]
            contig_acc[contig_id] = acc1

        acc2 = genome_acc.get(genome_id)
        if acc2:
            acc2[og_genome_depth] += depth
            acc2[og_genome_covered_bases] += 1
        else:
            acc2 = [depth, 1]
            genome_acc[genome_id] = acc2

        #if line < 10:
        #    tsprint(("\n" + json.dumps(dict(zip(columns.keys(), row)), indent=4)).replace("\n", "\n" + sample_pileup_path + ": "))

        A, C, G, T = row[c_A], row[c_C], row[c_G], row[c_T]
        number_alleles = 0
        nonzero_allele_index = 4
        for i, nt_count in enumerate((A, C, G, T)):
            if nt_count / depth >= param.MIN_ALLELE_FREQUENCY_WITHIN_SAMPLE:
                number_alleles += 1
                nonzero_allele_index = i
        # Doesn't matter which of the at-most-2 nonzero alleles we record here,
        # because their frequencies add up to 1 and we can always infer the
        # letter of the other one later, so one can be derived from the other.
        nonzero_allele = "ACGTN"[nonzero_allele_index]
        nonzero_allele_count = row[columns[nonzero_allele]]

        # Within-sample filter, partition the pileup for the threads
        if number_alleles > 2:
            continue
        if depth < param.MIN_DEPTH:
            continue

        assert len(row) == c_site_id
        row.append(site_id)
        assert len(row) == c_genome_id
        row.append(genome_id)
        assert len(row) == c_number_alleles
        row.append(number_alleles)
        assert len(row) == c_nz_allele
        row.append(nonzero_allele)
        assert len(row) == c_nz_allele_count
        row.append(nonzero_allele_count)

        sites[site_id] = row

    # print_top(contig_depth)
    # print_top(contig_covered_bases)
    # print_top(genome_covered_bases)

    with open(banded_output_path, "w") as o1:
        o1.write("\t".join(["site_id", "depth", "A", "C", "G", "T", "nz_allele", "nz_allele_count"]) + "\n")
        output_sites = 0
        for site_id, row in sites.items():
            if genome_accumulator[sample_name][row[c_genome_id]][og_genome_covered_bases] < param.MIN_GENOME_COVERED_BASES:
                continue
            o1.write("\t".join([row[c_site_id], str(row[c_depth]), str(row[c_A]), str(row[c_C]), str(row[c_G]), str(row[c_T]), row[c_nz_allele], str(row[c_nz_allele_count])]) + "\n")
            output_sites += 1

    t_end = time.time()
    #tsprint(f"{sample_name}_tid{thread_id}: Output {output_sites} sites passing filters, out of {len(sites)} total sites.  Pass rate: {output_sites/len(sites)*100:3.1f} percent.")
    tsprint(f"{sample_name}_tid{thread_id}: Run time {t_end - t_start} seconds, or {sites_count/(t_end - t_start):.1f} sites per second.")
    return "it worked"


def process_worker(args):
    sample_list_file, sample_file_names, num_threads, thread_id = args
    t_start = time.time()

    contig_accumulator = defaultdict(dict)
    genome_accumulator = defaultdict(dict)

    for sample_index, sample_pileup_path in enumerate(sample_file_names):
        process(sample_pileup_path, num_threads, thread_id, contig_accumulator, genome_accumulator)

    output_path_contig_stats = f"banded/band{thread_id}.contig_stats.tsv"
    output_path_genome_stats = f"banded/band{thread_id}.genome_stats.tsv"

    with open(output_path_contig_stats, "w") as o2:
        o2.write("sample_name\tcontig_id\tgenome_id\tcontig_total_depth\tcontig_covered_bases\n")
        for sample_name, contig_acc in contig_accumulator.items():
            for contig_id, contig_info in contig_acc.items():
                contig_total_depth, contig_covered_bases = contig_info
                genome_id = contig_id.split("_", 1)[0]
                o2.write(f"{sample_name}\t{contig_id}\t{genome_id}\t{contig_total_depth}\t{contig_covered_bases}\n")

    with open(output_path_genome_stats, "w") as o3:
        o3.write("sample_name\tgenome_id\tgenome_total_depth\tgenome_covered_bases\n")
        for sample_name, genome_acc in genome_accumulator.items():
            for genome_id, genome_info in genome_acc.items():
                genome_total_depth, genome_covered_bases = genome_info
                o3.write(f"{sample_name}\t{genome_id}\t{genome_total_depth}\t{genome_covered_bases}\n")

    t_end = time.time()
    tsprint(f"THREAD {thread_id}: Run time {t_end - t_start} seconds.")
    return "it worked"


def main():
    assert len(sys.argv) > 1
    sample_list_file = sys.argv[1]
    contig_accumulator = defaultdict(dict)
    genome_accumulator = defaultdict(dict)

    with open(sample_list_file, "r") as slf:
        sample_file_names = [line.strip() for line in slf]

    t_start = time.time()
    mp = multiprocessing.Pool(param.THREADS)
    results = mp.map(process_worker, [(sample_list_file, sample_file_names, param.THREADS, thread_id) for thread_id in range(param.THREADS)])
    t_end = time.time()

    tsprint(f"ALL THREADS:  Run time {t_end - t_start} seconds.")
    assert all(s == "it worked" for s in results)


if __name__ == "__main__":
    main()
