#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
import os
from collections import defaultdict

import param
from util import *


def accumulate(accumulator, sample_file_names, sample_brief_names, sample_index, thread_id):
    sample_pileup_path = sample_file_names[sample_index]
    sample_name = sample_brief_names[sample_index]
    samples_count = len(sample_file_names)
    band = sample_pileup_path.rsplit(".", 2)[1]
    paramstr = f"dp{param.MIN_DEPTH}.gcb{param.MIN_GENOME_COVERED_BASES}.{band}"

    input_path_contig_stats = f"banded/{band}.contig_stats.tsv"
    input_path_genome_stats = f"banded/{band}.genome_stats.tsv"

    # Load genome stats
    table_iterator = parse_table(tsv_rows(input_path_genome_stats), param.schema_genome_stats)
    columns = next(table_iterator)
    gs_sample_name = columns["sample_name"]
    gs_genome_id = columns["genome_id"]
    gs_total_depth = columns["total_depth"]
    gs_covered_bases = columns["covered_bases"]
    gs_coverage = columns["coverage"] = len(columns)
    genome_stats = defaultdict(dict)
    for line, row in enumerate(table_iterator):
        row.append(row[gs_total_depth] / row[gs_covered_bases])
        sname = row[gs_sample_name]
        genome_id = row[gs_genome_id]
        genome_stats[sname][genome_id] = row

    # Load contig stats
    table_iterator = parse_table(tsv_rows(input_path_contig_stats), param.schema_contig_stats)
    columns = next(table_iterator)
    cs_sample_name = columns["sample_name"]
    cs_contig_id = columns["contig_id"]
    cs_total_depth = columns["total_depth"]
    cs_covered_bases = columns["covered_bases"]
    cs_coverage = columns["coverage"] = len(columns)
    contig_stats = defaultdict(dict)
    for line, row in enumerate(table_iterator):
        row.append(row[cs_total_depth] / row[cs_covered_bases])
        sname = row[cs_sample_name]
        contig_id = row[cs_contig_id]
        contig_stats[sname][contig_id] = row

    # Read banded pileup files
    #table_iterator = parse_table(tsv_rows_slice2(sample_pileup_path, num_threads, thread_id), param.sample_pileup_schema_banded)
    table_iterator = parse_table(tsv_rows_slice_contig(sample_pileup_path, thread_id), param.sample_pileup_schema_banded_v2)
    columns = next(table_iterator)

    # Get integer keys for columns
    c_genome_id = columns["genome_id"]
    c_ref_id = columns["ref_id"]
    c_ref_pos = columns["ref_pos"]
    c_ref_allele = columns["ref_allele"]
    c_depth = columns["depth"]
    c_A = columns["A"]
    c_C = columns["C"]
    c_G = columns["G"]
    c_T = columns["T"]
    c_nz_allele = columns["nz_allele"]
    c_nz_allele_count = columns["nz_allele_count"]

    # Output column indices
    s_A, s_C, s_G, s_T, s_sample_count, s_scA, s_scC, s_scG, s_scT = range(9)

    for line, row in enumerate(table_iterator):

        if line % (1000*1000) == 0:
            tsprint(f"{sample_pileup_path}_THREAD{thread_id}_SAMPLEIDX{sample_index}: Processing {line}.")
        if line == param.MAX_LINES:
            break

        # Unpack frequently accessed columns.
        genome_id = row[c_genome_id]
        contig_id = row[c_ref_id]
        ref_pos = row[c_ref_pos]
        ref_allele = row[c_ref_allele]
        depth = row[c_depth]
        nz_allele = row[c_nz_allele]
        nz_allele_count = row[c_nz_allele_count]
        A, C, G, T = row[c_A], row[c_C], row[c_G], row[c_T]

        # Compute derived columns.
        site_id = f"{contig_id}|{ref_pos}|{ref_allele}"
        nz_allele_freq = nz_allele_count / depth
        site_ratio = depth / contig_stats[sample_name][contig_id][cs_coverage]
        genome_coverage = genome_stats[sample_name][genome_id][gs_coverage]

        # Filter.
        if genome_coverage < param.MIN_GENOME_COVERAGE:
            continue
        if depth < param.MIN_DEPTH_SNP:
            continue
        if site_ratio > param.MAX_SITE_RATIO:
            continue

        # Sample count for ACGT
        sc_ACGT = [0, 0, 0, 0]
        for i, nt_count in enumerate((A, C, G, T)):
            if nt_count / depth >= param.MIN_ALLELE_FREQUENCY_WITHIN_SAMPLE:
                sc_ACGT[i] = 1

        # Aggregate.
        genome_acc = accumulator[genome_id]
        acc = genome_acc.get(site_id)
        if acc:
            acc[s_A] += A
            acc[s_C] += C
            acc[s_G] += G
            acc[s_T] += T
            acc[s_sample_count] += 1
            acc[s_scA] += sc_ACGT[0]
            acc[s_scC] += sc_ACGT[1]
            acc[s_scG] += sc_ACGT[2]
            acc[s_scT] += sc_ACGT[3]
        else:
            acc = [A, C, G, T, 1, sc_ACGT[0], sc_ACGT[1], sc_ACGT[2], sc_ACGT[3]] + ([('N', 0)] * samples_count)
            genome_acc[site_id] = acc
        # This isn't being accumulated across samples;  we are just remembering the value from each sample.
        assert acc[9 + sample_index] == ('N', 0) and nz_allele != 'N'
        acc[9 + sample_index] = (nz_allele, nz_allele_freq)


def filter2(accumulator, sample_list_file, sample_brief_names):

    outpref = os.path.basename(sample_list_file).rsplit(".", 1)[0]
    for genome_id, genome_acc in accumulator.items():

        output_sites = f"banded/accumulators_{outpref}.gid_{genome_id}.dp_{param.MIN_DEPTH_SNP}.mgc_{param.MIN_GENOME_COVERAGE}.sr_{param.MAX_SITE_RATIO}.tsv"

        with open(output_sites, "w") as out_sites:
            out_sites.write("site_id\tA\tC\tG\tT\tsample_count\tscA\tscC\tscG\tscT\n")
            out_sites.write("\t".join(["major_allele", "minor_allele"] + sample_brief_names) + "\n")
            for site_id, site_info in genome_acc.items():
                A, C, G, T, sample_count, scA, scC, scG, scT = site_info[:9]
                depth = A + C + G + T
                all_alleles = ((A, 'A'), (C, 'C'), (G, 'G'), (T, 'T'))
                alleles_above_cutoff = tuple(al for al in all_alleles if al[0] / depth >= param.MIN_ALLELE_FREQUECY_ACROSS_SAMPLES)
                # Keep only bi-allelic and mono-allelic sites.
                if 1 <= len(alleles_above_cutoff) <= 2:
                    # In the event of a tie -- biallelic site with 50/50 freq split -- the allele declared major is
                    # the one that comes later in the "ACGT" lexicographic order.
                    alleles_above_cutoff = sorted(alleles_above_cutoff, reverse=True)
                    major_allele = alleles_above_cutoff[0][1]
                    minor_allele = alleles_above_cutoff[-1][1]  # for mono-allelic sites, same as major allele
                    out_sites.write(f"{site_id}\t{A}\t{C}\t{G}\t{T}\t{sample_count}\t{scA}\t{scC}\t{scG}\t{scT}\t")
                    major_allele_freqs_by_sample = "\t".join(
                        "{:.3f}".format(-1.0 if allele == 'N' else (freq if allele==major_allele else 1.0 - freq))
                        for allele, freq in site_info[9:])
                    out_sites.write(major_allele + "\t" + minor_allele + "\t" + major_allele_freqs_by_sample + "\n")


def process_worker(args):
    sample_list_file, sample_file_names, thread_id = args
    t_start = time.time()
    accumulator = defaultdict(dict)
    sample_brief_names = [os.path.basename(sfn).split(".", 1)[0] for sfn in sample_file_names]

    for sample_index, sample_pileup_path in enumerate(sample_file_names):
        accumulate(accumulator, sample_file_names, sample_brief_names, sample_index, thread_id)
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

    # Handle the multi-process here
    num_threads = param.THREADS
    threads_banded = defaultdict(dict)
    for thread_id in range(num_threads):
        genomes = [g for g in list(set(param.CONTIGS.values())) if int(g[4:6]) % num_threads == thread_id]
        contigs = [c for c,g in param.CONTIGS.items() if g in genomes]
        threads_banded[thread_id]['genomes'] = list(genomes)
        threads_banded[thread_id]['contigs'] = list(contigs)
        output_path_contig_list = f"banded/band{thread_id}.contig_lists.txt"
        with open(output_path_contig_list, 'w') as ot:
            for contig in contigs:
                ot.write(f"{contig}\n")

    t_start = time.time()
    mp = multiprocessing.Pool(param.THREADS)
    results = mp.map(process_worker, [(sample_list_file, sample_file_names, thread_id) for thread_id in range(param.THREADS)])
    t_end = time.time()

    tsprint(f"ALL THREADS:  Run time {t_end - t_start} seconds.")
    assert all(s == "it worked" for s in results)


if __name__ == "__main__":
    main()
