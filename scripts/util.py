#!/usr/bin/env python3
import sys
import json
import multiprocessing
import time


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


def tsv_rows_slice(path, num_threads, thread_id):
    # TODO:  Support s3 and compressed files.
    with open(path, "r") as stream:
        yield next(stream).rstrip("\n").split("\t")
        for line in stream:
            if int(line[4:6]) % num_threads == thread_id:
                yield line.rstrip("\n").split("\t")


def print_top(counters, how_many=5):
    print(json.dumps(sorted(((depth, contig_id) for contig_id, depth in counters.items()), reverse=True)[:how_many], indent=4))


def chomp(s, ending):
    assert s.endswith(ending)
    return s[:-len(ending)]


if __name__ == "__main__":
    tsprint("Hello from xsnp utils.")
