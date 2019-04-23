#!/usr/bin/env python
import json
from collections import defaultdict

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
    for values in rows:
        assert len(headers) == len(values)
        yield dict(zip(headers, [f(v) for f,v in zip(functions, values)]))

def tsv_rows(path):
    # TODO:  Support s3 and compressed files.
    with open(path, "r") as stream:
        for line in stream:
            yield line.rstrip("\n").split("\t")

if __name__ == "__main__":
    site_depth = defaultdict(int)
    contig_depth = defaultdict(int)
    covered_bases  = defaultdict(int)
    for line, row in enumerate(parse_table(tsv_rows("SRS023176.pileup"), schema)):
        # if line == 100000:
        #    break
        contig_id = row["ref_id"]
        row_depth = row["depth"]
        site_id = contig_id + "|" + str(row["ref_pos"]) + "|" + str(row["ref_allele"])
        site_depth[site_id] += row_depth
        contig_depth[contig_id] += row_depth
        covered_bases[contig_id] += 1
        if line < 10:
            print(json.dumps(row, indent=4))
    # print(json.dumps(sorted(((depth, site_id) for site_id, depth in site_depth.items()), reverse=True)[:10], indent=4))
    print(json.dumps(sorted(((depth, contig_id) for contig_id, depth in contig_depth.items()), reverse=True)[:10], indent=4))
    print(json.dumps(sorted(((site_count, contig_id) for contig_id, site_count in covered_bases.items()), reverse=True)[:10], indent=4))
