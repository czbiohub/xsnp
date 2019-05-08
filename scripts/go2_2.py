#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
from collections import defaultdict

import param
from util import *


def accumulate(accumulator, sample_file_names, sample_index, num_threads, thread_id, contig_stats, genome_stats):

    sample_pileup_path = sample_file_names[sample_index]
    sample_name = chomp(sample_pileup_path, ".pileup")

    paramstr = f"gcb{param.MIN_GENOME_COVERED_BASES}.dp{param.MIN_DEPTH}.tid{thread_id}"
    sample_pileup_tid_path = f"intermediate/{sample_name}.sites.{paramstr}.tsv"
    tsprint(f"{sample_pileup_tid_path}: Processing sample pileup path")

    table_iterator = parse_table(tsv_rows(sample_pileup_tid_path), param.sample_pileup_id_schema)
    #table_iterator = parse_table(tsv_rows_slice(sample_pileup_path, num_threads, thread_id), param.sample_pileup_schema)
    columns = next(table_iterator)

    # Get integer keys for columns
    c_site_id = columns["site_id"]
    c_depth = columns["depth"]
    c_A = columns["A"]
    c_C = columns["C"]
    c_G = columns["G"]
    c_T = columns["T"]

    # output column indices
    s_A, s_C, s_G, s_T, s_sample_count = range(5)

    for line, row in enumerate(table_iterator):

        if line % (1000*1000 == 0:
            tsprint(f"{sample_pileup_path}:{thread_id}:{sample_index} Processing {line}.")
        if line == param.MAX_LINES:
            break

        # Unpack frequently accessed columns.
        site_id = row[c_site_id]
        depth = row[c_depth]
        A, C, G, T = row[c_A], row[c_C], row[c_G], row[c_T]

        # Compute derived columns.
        genome_id = site_id.split("_", 1)[0]
        contig_id = site_id.split("|", 1)[0]

        site_ratio = depth / contig_stats[sample_name][contig_id][-1]
        genome_coverage = genome_stats[sample_name][genome_id][-1]

        # Filter.
        if genome_coverage < param.MIN_GENOME_COVERAGE:
            continue
        if site_ratio > param.MAX_SITE_RATIO:
             continue

        # Aggregate.
        genome_acc = accumulator[genome_id]
        acc = genome_acc.get(site_id)
        if acc:
            acc[s_A] += A
            acc[s_C] += C
            acc[s_G] += G
            acc[s_T] += T
            acc[s_sample_count] += 1
        else:
            acc = [A, C, G, T, 1]
            genome_acc[site_id] = acc


def filter2(accumulator, sample_list_file, sample_brief_names):

    outpref = sample_list_file.rsplit(".", 1)[0]
    for genome_id, genome_acc in accumulator.items():

        output_sites = f"accumulators_{outpref}.gid_{genome_id}.sr_{param.MAX_SITE_RATIO}.mgc_{param.MIN_GENOME_COVERAGE}.tsv"

        if len(genome_acc) > 0:
            with open(output_sites, "w") as out_sites:
                out_sites.write("site_id\tA\tC\tG\tT\tsample_count\tmajor_allele\tminor_allele\t")
                ## read them and filter them
                out_sites.write("\t".join(sample_brief_names) + "\n")
                for site_id, site_info in genome_acc.items():
                    A, C, G, T, sample_count = site_info[:5]
                    depth = A + C + G + T
                    all_alleles = ((A, 'A'), (C, 'C'), (G, 'G'), (T, 'T'))
                    alleles_above_cutoff = tuple(al for al in all_alleles if al[0] / depth >= param.MIN_ALLELE_FREQUECY_ACROSS_SAMPLES)
                    # Keep only bi-allelic and mono-allelic sites across samples.
                    if 1 <= len(alleles_above_cutoff) <= 2:
                        # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
                        # the one that comes later in the "ACGT" lexicographic order.
                        alleles_above_cutoff = sorted(alleles_above_cutoff, reverse=True)
                        major_allele = alleles_above_cutoff[0][1]
                        minor_allele = alleles_above_cutoff[-1][1]  # for mono-allelic sites, same as major allele
                        out_sites.write(f"{site_id}\t{A}\t{C}\t{G}\t{T}\t{sample_count}\t{major_allele}\t{minor_allele}\t")
                        major_allele_freqs_by_sample = "\t".join(
                            "{:.3f}".format(0.0 if allele == 'N' else (freq if allele==major_allele else 1.0 - freq))
                            for allele, freq in site_info[5:])
                        out_sites.write(major_allele_freqs_by_sample + "\n")


def process_worker(args):
    sample_list_file, sample_file_names, num_threads, thread_id = args
    t_start = time.time()


    accumulator = defaultdict(dict)
    sample_brief_names = [chomp(sfn, ".pileup") for sfn in sample_file_names]
    for sample_index, sample_pileup_path in enumerate(sample_file_names):
        accumulate(accumulator, sample_file_names, sample_index, num_threads, thread_id, contig_stats, genome_stats)

    # Read and filter
    samples_count = len(sample_brief_names)
    for sample_index, sample_name in enumerate(sample_brief_names):
        paramstr = f"gcb{param.MIN_GENOME_COVERED_BASES}.dp{param.MIN_DEPTH}.tid{thread_id}"
        sample_pileup_tid_path = f"intermediate/{sample_name}.sites.{paramstr}.tsv"

        table_iterator = parse_table(tsv_rows(sample_pileup_tid_path), param.sample_pileup_id_schema)
        columns = next(table_iterator)
        # Get integer keys for columns
        c_site_id = columns["site_id"]
        c_nz_allele = columns["nz_allele"]
        c_nz_allele_freq = columns["nz_allele_freq"]
        for line, row in enumerate(table_iterator):
            # Unpack frequently accessed columns.
            site_id = row[c_site_id]
            nz_allele = row[c_nz_allele]
            nz_allele_freq = row[c_nz_allele_freq]
            # Compute derived columns.
            genome_id = site_id.split("_", 1)[0]
            genome_acc = accumulator[genome_id]
            acc = genome_acc.get(site_id)
            if acc:
                if len(acc) == 5:
                    acc = acc + ([('N', 0)] * samples_count)
                    genome_acc[site_id] = acc
                # This isn't being accumulated across samples;  we are just remembering the value from each sample for each site.
                assert acc[5 + sample_index] == ('N', 0) and nz_allele != 'N'
                acc[5 + sample_index] = (nz_allele, nz_allele_freq)

    # Filter
    filter2(accumulator, sample_list_file, sample_brief_names)
    t_end = time.time()
    tsprint(f"THREAD {thread_id}: Run time {t_end - t_start} seconds.")
    return "it worked"


def main():
    assert len(sys.argv) > 1
    accumulator = defaultdict(dict)
    sample_list_file = sys.argv[1]
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
