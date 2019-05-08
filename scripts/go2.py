#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time
from collections import defaultdict

import param
from util import *


def accumulate(accumulator, sample_file_names, sample_brief_names, sample_index, num_threads, thread_id):

    sample_pileup_path = sample_file_names[sample_index]
    #tsprint(f"{sample_pileup_path}: Processing sample pileup path")
    sample_name = sample_brief_names[sample_index]
    samples_count = len(sample_file_names)
    band = sample_pileup_path.rsplit(".", 2)[1]
    paramstr = f"gcb{param.MIN_GENOME_COVERED_BASES}.dp{param.MIN_DEPTH}.{band}"

    input_path_contig_stats = f"{band}.contig_stats.tsv"
    input_path_genome_stats = f"{band}.genome_stats.tsv"

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

    # Read banhded pileup files
    table_iterator = parse_table(tsv_rows_slice2(sample_pileup_path, num_threads, thread_id), param.sample_pileup_schema_banded)
    columns = next(table_iterator)

    # Get integer keys for columns
    c_site_id = columns["site_id"]
    c_depth = columns["depth"]
    c_A = columns["A"]
    c_C = columns["C"]
    c_G = columns["G"]
    c_T = columns["T"]
    c_nz_allele = columns["nz_allele"]
    c_nz_allele_count = columns["nz_allele_count"]

    # Output column indices
    s_A, s_C, s_G, s_T, s_sample_count = range(5)

    for line, row in enumerate(table_iterator):

        if line % (1000*1000) == 0:
            tsprint(f"{sample_pileup_path}:tid{thread_id}:sidx{sample_index} Processing {line}.")
        if line == param.MAX_LINES:
            break

        # Unpack frequently accessed columns.
        site_id = row[c_site_id]
        depth = row[c_depth]
        nz_allele = row[c_nz_allele]
        nz_allele_count = row[c_nz_allele_count]
        A, C, G, T = row[c_A], row[c_C], row[c_G], row[c_T]

        # Compute derived columns.
        genome_id = site_id.split("_", 1)[0]
        contig_id = site_id.split("|", 1)[0]
        nz_allele_freq = nz_allele_count / depth

        site_ratio = depth / contig_stats[sample_name][contig_id][cs_coverage]
        genome_coverage = genome_stats[sample_name][genome_id][gs_coverage]

        # Filter.
        if depth < param.MIN_DEPTH_SNP:
            continue
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
            acc = [A, C, G, T, 1] + ([('N', 0)] * samples_count)
            genome_acc[site_id] = acc
        # This isn't being accumulated across samples;  we are just remembering the value from each sample.
        assert acc[5 + sample_index] == ('N', 0) and nz_allele != 'N'
        acc[5 + sample_index] = (nz_allele, nz_allele_freq)


def filter2(accumulator, sample_list_file, sample_brief_names):

    outpref = sample_list_file.rsplit(".", 1)[0]
    for genome_id, genome_acc in accumulator.items():

        output_sites = f"accumulators_{outpref}.gid_{genome_id}.sr_{param.MAX_SITE_RATIO}.dp_{param.MIN_DEPTH_SNP}.mgc_{param.MIN_GENOME_COVERAGE}.tsv"

        with open(output_sites, "w") as out_sites:
            out_sites.write("site_id\tA\tC\tG\tT\tsample_count\t")
            out_sites.write("\t".join(["major_allele", "minor_allele"] + sample_brief_names) + "\n")
            for site_id, site_info in genome_acc.items():
                A, C, G, T, sample_count = site_info[:5]
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
                    out_sites.write(f"{site_id}\t{A}\t{C}\t{G}\t{T}\t{sample_count}\t")
                    major_allele_freqs_by_sample = "\t".join(
                        "{:.3f}".format(0.0 if allele == 'N' else (freq if allele==major_allele else 1.0 - freq))
                        for allele, freq in site_info[5:])
                    out_sites.write(major_allele + "\t" + minor_allele + "\t" + major_allele_freqs_by_sample + "\n")


def process_worker(args):
    sample_list_file, sample_file_names, num_threads, thread_id = args
    t_start = time.time()
    accumulator = defaultdict(dict)
    sample_brief_names = [sfn.split(".", 1)[0] for sfn in sample_file_names]
    for sample_index, sample_pileup_path in enumerate(sample_file_names):
        accumulate(accumulator, sample_file_names, sample_brief_names, sample_index, num_threads, thread_id)
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
