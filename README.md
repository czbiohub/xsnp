# xsnp
Cross-Sample SNP Analysis Tools

Broadly, the input consists of aignment pileup data for a large number of samples and species.   The output is a compact representation of the input that can be loaded into R to enable downstream analysis with PCA, tSNE, clusterig, tree building, etc.

# 1.  Input:  Table `pileup` with schema
```
    sample_id       "SRS050026"
    contig_id       "261195|accn|NZ_ADKO01000092"
    genome_id       261195
    ref_pos         396
    ref_allele      "C"                     *-------------------------------------------*
    depth           3450                    *    KEY:  (sample_id, contig_id, ref_pos)  *
    count_a         3                       *-------------------------------------------*
    count_c         3446
    count_g         1
    count_t         0
```

# 2.  Input:  Interpetation

A reference genome is comprised of one or more contigs, where
```
contig_id.starts_with(genome_id)
```

The pair `(contig_id, ref_pos)` constitute a genomic `site_id`.  The nucleotide at that site in the reference genome is `ref_allele`.  It need not be `"A"`, `"C"`, `"G"`, or `"T"`; might be `"N"`, `"Y"`, `"R"`, etc.
```
0 <= ref_pos < len(ref_contig)

ref_allele = ref_contig[ref_pos]
```

Pileup `depth` is the number of reads from `sample_id` aligned over `site_id`, including only reads that have `"A"`, `"C"`, `"G"`, or `"T"` at that site (even if the reference allele is non-ACGT).  Among those reads, `count_a` have `"A"` at the site, `count_c` have `"C"`, etc.
```
depth = count_a + count_c + count_g + count_t

depth > 0
```
These counts are a function of `(sample_id, contig_id, ref_pos)`, the key to the pileup schema.


# 3.  Notation


Let `t.k` denote the row in table `t` that has key valye `k`.

Let `t.k.c` denote the value in column `c` of row `t.k`.

Let us omit `t` when it is clear from context.  This is natural when referring to object properties, for example
```
    genome_id.length
        
    (sample_id, site_id).ref_allele
```
instead of the needlessly verbose
```
    ref_genomes.genome_id.length
    
    pileup.(sample_id, site_id).ref_allele    
```
The reference allele depends only on `site_id`, not `sample_id`.  It may be abbreviated to
```
    site_id.ref_allele
```
A more disciplinned approach would pull out data that only depends on the site into a separate `ref_genomes` table.  In practice we have so little of this data (literally just the `ref_allele` and `genome_id` columns) that it isn't worth pulling into its own table for the purpose of this particular computation.

# 4.  The `matches` family of predicates

The following predicates are defined for any row that has the necessary columns, regardless of the table that row is from.

```
    row.matches(genome_id, sample_id) := (row.genome_id == genome_id) AND (row.sample_id == sample_id)
    
    row.matches(contig_id, sample_id) := (row.contig_id == contig_id) AND (row.sample_id == sample_id)
    
    row.matches(contig_id, ref_pos) := (row.contig_id == contig_id) AND (row.ref_pos == ref_pos)
```


# 5.  Derived table `contig_stats`
Table `contig_stats` is computed in the first pass over the pileup data, and has schema
```
    genome_id
    contig_id
    sample_id
    depth
    covered_bases
```
where
```
    depth := SUM(pileup_row.depth) WHERE pileup_row.matches(contig_id, sample_id)

    covered_bases := SUM(1) WHERE pileup_row.matches(contig_id, sample_id)

    coverage := depth / covered_bases
```

# 6.  Derived table `genome_stats`
Table `genome_stats` is also computed in the first pass, along with `contig_stats`.  It is defined very similarly to `contig_stats`, with schema
```
    genome_id
    sample_id
    depth
    covered_bases
```
where
```
    depth := SUM(contig_stats_row.depth) WHERE contig_stats_row.matches(genome_id, sample_id)

    covered_bases := SUM(1) WHERE contig_stats_row.matches(genome_id, sample_id)
```

# 7. Derived columns added to table `pileup`
A second pass over the input data extends table `pileup` with the following derived columns:
```
    WITHIN pileup_row:

        num_alleles := SUM(1) FOR c in [count_a, count_c, count_g, count_t] WHERE c / depth > MIN_ALLELE_FREQUENCY_WITHIN_SAMPLE

        of_interest := num_alleles in (1, 2)
                       AND depth >= MIN_DEPTH
                       AND depth <= MAX_SITE_RATIO * contig_stats.(contig_id, sample_id).coverage
                       AND genome_stats.(genome_id, sample_id).coverage >= MIN_GENOME_COVERAGE
                       AND genome_stats.(genome_id, sample_id).covered_bases >= MIN_GENOME_COVERED_BASES
```
An alternative representation would be to create a separate new table with the same key as table `pileup`.

# 8.  Derived table `prelim_results`
The second pass also creates a new table with schema
```
   contig_id
   ref_pos
   A
   C
   G
   T
   sample_count
   major_allele                   *-----------------------------*
   minor_allele                   *  KEY: (contig_id, ref_pos)  *
   depth                          *-----------------------------*
   num_alleles
   of_interest
```
and columns defined as follows
```
    A, C, G, T := SUM(pileup_row.(count_a, count_c, count_g, count_t)) WHERE pileup_row.of_interest AND pileup_row.matches(contig_id, ref_pos)
   
    depth := A + C + G + T
    
    sample_count := COUNT DISTINCT pileup_row.sample_id WHERE pileup_row.of_interest AND pileup_row.matches(contig_id, ref_pos)
   
    num_alleles := SUM(1) FOR c in [A, C, G, T] if c / depth >= MIN_ALLELE_FREQUECY_ACROSS_SAMPLES
   
    of_interest := num_alleles in (1, 2)
   
    _alleles_ := sorted_desecending((A, 'A'), (C, 'C'), (G, 'G'), (T, 'T))
   
    major_allele := _alleles_[0][1]
   
    minor_allele := _alleles_[1 % num_alleles][1]
```

# 9.  Cross-sample counts
We are interested in a subset of the prelim_results table.
```
    SELECT * FROM prelim_results WHERE of_interest
```


# 10.  Within-sample counts
Lastly, we add two more derived columns to table `pileup`, defined only where
```prelim_results.(contig_id, ref_pos).of_interest```
These are
```
    WITHIN pileup_row:
    
        major_allele := prelim_results.(contig_id, ref_pos).major_allele
    
        major_allele_frequency := count_{major_allele} / depth
```
This can be a byproduct of pass 2, or a separate pass 3, depending on implementation and available resources.
