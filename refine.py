#!/usr/bin/env python3
import argparse
import math

COMMENT_LINE='#'

def refine(file, skip=0):
    agg_mean = 0
    agg_err = 0
    agg_w = 0

    for line in file:
        line = line.strip()
        if line.startswith(COMMENT_LINE):
            continue
        if skip > 0:
            skip -= 1
            continue

        mean, err, w = line.split(',')
        mean = float(mean)
        err = float(err)
        w = int(w)

        agg_mean += w * mean
        agg_err += w * w * err * err
        agg_w += w

    agg_mean /= agg_w
    agg_err = math.sqrt(agg_err)
    agg_err /= agg_w

    inv_mean = 1.0 / agg_mean
    inv_err = agg_err * inv_mean * inv_mean

    return inv_mean, inv_err

def csv_text(x):
    allowed = set(list("abcdefghijklmnopqrstuvwxyz0123456789-_ "))
    chars = set(list(x.lower()))
    if chars.issubset(allowed):
        return x
    return '"%s"' % x.encode('unicode_escape').replace(b'"', b'\\"').decode('ascii')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-s', '--skip', action='store', default=0, type=int,
                        help='number of output datums to skip')
    parser.add_argument('-c', '--csv', action='store_true', default=False)
    parser.add_argument('file', type=argparse.FileType('r'), nargs='*',
                        help='(clean) output file')

    args = parser.parse_args()
    skip = args.skip
    csv = args.csv

    for fa in args.file:
        with fa as file:
            try:
                inv_mean, inv_err = refine(file, skip)
                if csv: print("%s,ok,%r,%r"       % (csv_text(file.name), inv_mean, inv_err))
                else:   print("%s mean: %r (±%r)" % (file.name, inv_mean, inv_err))
            except Exception as e:
                if csv: print("%s,error,%s," % (csv_text(file.name), csv_text(repr(e))))
                else:   print("%s error: %r" % (file.name, e))