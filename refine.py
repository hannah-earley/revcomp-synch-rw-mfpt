#!/usr/bin/env python3
import argparse
import math

COMMENT_LINE='#'

class Experiment:
    def __init__(self):
        self.agg_mean = 0
        self.agg_err = 0
        self.agg_w = 0

    def insert(self, mean, err, w):
        mean = float(mean)
        err = float(err)
        w = int(w)

        self.agg_mean += w * mean
        self.agg_err += w * w * err * err
        self.agg_w += w

    def update(self, lines, skip=0):
        for line in lines:
            line = line.strip()
            if line.startswith(COMMENT_LINE):
                continue
            if skip > 0:
                skip -= 1
                continue
            self.insert(*line.split(','))

    def summarise(self):
        w = self.agg_w
        mean = self.agg_mean / w
        err = math.sqrt(self.agg_err) / w

        inv_mean = 1.0 / mean
        inv_err = err * inv_mean * inv_mean
        return inv_mean, inv_err

def refine(file, skip=0):
    expt = Experiment()
    expt.update(file, skip)
    return expt.summarise()

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