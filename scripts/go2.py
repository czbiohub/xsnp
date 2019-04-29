#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
from collections import defaultdict

import param
from util import *


def accumulate(accumulator, sample_pileup_path, num_threads, thread_id):

    tsprint(f"{sample_pileup_path}: Processing sample pileup path")
    sample_name = chomp(sample_pileup_path, ".pileup")
    paramstr = f"gcb{param.MIN_GENOME_COVERED_BASES}.dp{param.MIN_DEPTH}.sr{param.MAX_SITE_RATIO}"
    input_path_site_id = f"{sample_name}.sites.{paramstr}.tsv"
    input_path_contig_stats = f"{sample_name}.contig_stats.tsv"
    input_path_genome_stats = f"{sample_name}.genome_stats.tsv"

    # Load genome stats
    table_iterator = parse_table(tsv_rows(input_path_genome_stats), param.schema_genome_stats)
    columns = next(table_iterator)
    gs_genome_id = columns["genome_id"]
    gs_total_depth = columns["total_depth"]
    gs_covered_bases = columns["covered_bases"]
    gs_coverage = columns["coverage"] = len(columns)
    genome_stats = {}
    for line, row in enumerate(table_iterator):
        row.append(row[gs_total_depth] / row[gs_covered_bases])
        genome_id = row[gs_genome_id]
        genome_stats[genome_id] = row

    # Load contig stats
    table_iterator = parse_table(tsv_rows(input_path_contig_stats), param.schema_contig_stats)
    columns = next(table_iterator)
    cs_contig_id = columns["contig_id"]
    cs_total_depth = columns["total_depth"]
    cs_covered_bases = columns["covered_bases"]
    cs_coverage = columns["coverage"] = len(columns)
    contig_stats = {}
    for line, row in enumerate(table_iterator):
        row.append(row[cs_total_depth] / row[cs_covered_bases])
        contig_id = row[cs_contig_id]
        contig_stats[contig_id] = row

    #table_iterator = parse_table(tsv_rows(sample_pileup_path), param.sample_pileup_schema)
    table_iterator = parse_table(tsv_rows_slice(sample_pileup_path, num_threads, thread_id), param.sample_pileup_schema)
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

    s_sample_count = 4

    for line, row in enumerate(table_iterator):

        if line % (1000*1000) == 0:
            tsprint(f"{sample_pileup_path}:{thread_id}: Processing {line}.")
        if line == param.MAX_LINES:
            break
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
            if row[nt_count] / row_depth >= param.MIN_ALLELE_FREQUENCY_WITHIN_SAMPLE:
                number_alleles += 1
        assert len(row) == c_number_alleles
        row.append(number_alleles)

        #if line < 10:
        #    tsprint(("\n" + json.dumps(dict(zip(columns.keys(), row)), indent=4)).replace("\n", "\n" + sample_pileup_path + ": "))

        if number_alleles > 2:
            continue

        if genome_stats[genome_id][gs_coverage] < param.MIN_GENOME_COVERAGE:
            continue

        if row[c_depth] < param.MIN_DEPTH:
            continue

        if genome_stats[genome_id][gs_covered_bases] < param.MIN_GENOME_COVERED_BASES:
            continue

        contig_id = row[c_ref_id]
        site_ratio = row[c_depth] / contig_stats[contig_id][cs_coverage]
        if site_ratio > param.MAX_SITE_RATIO:
             continue

        # Aggregate
        genome_acc = accumulator[genome_id]
        acc = genome_acc.get(site_id)
        if not acc:
            acc = [0, 0, 0, 0, 0]
            genome_acc[site_id] = acc

        acc[s_sample_count] += 1

        for i, nt_count in enumerate((c_A, c_C, c_G, c_T)):
            acc[i] += row[nt_count]


def filter2(accumulator, sample_list_file):

    outpref = sample_list_file.rsplit(".", 1)[0]
    for genome_id, genome_acc in accumulator.items():

        output_sites = f"accumulators_{outpref}.gid_{genome_id}.sr_{param.MAX_SITE_RATIO}.mgc_{param.MIN_GENOME_COVERAGE}.tsv"

        with open(output_sites, "w") as out_sites:
            out_sites.write("site_id\tA\tC\tG\tT\tsample_count\n")  # hi
            for site_id, site_info in genome_acc.items():
                A, C, G, T, sample_count = site_info
                depth = A + C + G + T
                num_alleles = 0
                for ac in (A, C, G, T):
                    if ac / depth >= param.MIN_GENOME_COVERAGE:
                        num_alleles += 1
                if num_alleles > 2:
                    continue
                out_sites.write(f"{site_id}\t{A}\t{C}\t{G}\t{T}\t{sample_count}\n")


def process_worker(args):
    sample_list_file, num_threads, thread_id = args
    t_start = time.time()
    accumulator = defaultdict(dict)
    with open(sample_list_file, "r") as slf:
        for sample_pileup in slf:
            accumulate(accumulator, sample_pileup.strip(), num_threads, thread_id)
    filter2(accumulator, sample_list_file)
    t_end = time.time()
    tsprint(f"THREAD {thread_id}: Run time {t_end - t_start} seconds.")
    return "it worked"


def main():
    assert len(sys.argv) > 1
    accumulator = defaultdict(dict)
    sample_list_file = sys.argv[1]
    t_start = time.time()
    mp = multiprocessing.Pool(param.THREADS)
    results = mp.map(process_worker, [(sample_list_file, param.THREADS, thread_id) for thread_id in range(param.THREADS)])
    t_end = time.time()
    tsprint(f"ALL THREADS:  Run time {t_end - t_start} seconds.")
    assert all(s == "it worked" for s in results)


if __name__ == "__main__":
    main()
