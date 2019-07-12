#!/usr/bin/env python3


THREADS = 24


# set to -1 to unlimit
MAX_LINES =  -1 #1 * 1000 * 100 #-1


MIN_DEPTH = 2
MIN_GENOME_COVERED_BASES = 10
MAX_SITE_RATIO = 5
MIN_ALLELE_FREQUENCY_WITHIN_SAMPLE = 0.01

MIN_DEPTH_SNP = 5
MIN_ALLELE_FREQUECY_ACROSS_SAMPLES = 0.05
MIN_GENOME_COVERAGE = 2.0

sample_pileup_schema = {
   "count_a": (int, "A"),
   "count_c": (int, "C"),
   "count_g": (int, "G"),
   "count_t": (int, "T"),
   "depth":   (int,),
   "ref_pos": (int,),
   "ref_id": (str, ),
   #"ref_id": (lambda ref_id: ref_id.replace("|", "_"),),
   "ref_allele": (str,)
}


schema_genome_stats = {
    "sample_name": (str,),
    "genome_id":(str,),
    "genome_total_depth": (int, "total_depth"),
    "genome_covered_bases": (int, "covered_bases")
}


schema_contig_stats = {
    "sample_name": (str,),
    "contig_id": (str,),
    "genome_id": (str,),
    "contig_total_depth": (int, "total_depth"),
    "contig_covered_bases": (int, "covered_bases")
}

sample_pileup_schema_banded_v1 = {
    "site_id": (str,),
    "depth": (int,),
    "A": (int,),
    "C": (int,),
    "G": (int,),
    "T": (int,),
    "nz_allele": (str,),
    "nz_allele_count": (int,)
}

sample_pileup_schema_banded_v2 = {
    "genome_id":(str,),
    "ref_id": (str, ),
    "ref_pos": (int,),
    "ref_allele": (str,),
    "depth": (int,),
    "A": (int,),
    "C": (int,),
    "G": (int,),
    "T": (int,),
    "number_alleles": (int,),
    "nz_allele": (str,),
    "nz_allele_count": (int,)
}

## 20190709: Add global dict mapping_contig_genome
CONTIGS = {}
with open("mapping_contig_genome.tsv") as f:
    for line in f:
        genome_id = line.rstrip('\n').split('\t')[0]
        contig_id = line.rstrip('\n').split('\t')[1]
        CONTIGS[contig_id] = str(genome_id)
