from __future__ import print_function
import math
import itertools as it
import array
from collections import defaultdict

import numpy as np
from scipy import interpolate
from matplotlib import ticker, pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

def maxend(fbed):
    emax = 0
    for toks in (x.rstrip().split("\t") for x in open(fbed)):
        try:
            emax = max(int(toks[2]), emax)
        except ValueError:
            continue
    return emax

def main(fbed, ped=None, prefix=None):
    e = maxend(fbed)
    d = {'1': np.zeros(e, dtype=np.uint16),
         '2': np.zeros(e, dtype=np.uint16)}

    lens, last_chrom = [], None

    sample_counts = defaultdict(lambda: defaultdict(int))
    sex = {x[1]: x[4] for x in (l.split("\t", 5) for l in open(ped))}

    for i, toks in enumerate(x.rstrip().split("\t") for x in open(fbed)):
        if i == 0:
            header = toks
            parent_id = toks.index('parent_id')
            continue
        pid = toks[parent_id]
        sample_counts[sex[pid]][pid] += 1

        if last_chrom != toks[0]:
            if last_chrom is not None:
                write(last_chrom, d, lens, prefix)
                d['1'][:] = 0
                d['2'][:] = 0
            lens = []
            last_chrom = toks[0]

        s, e = int(toks[1]), int(toks[2])
        lens.append(e - s)
        d[sex[pid]][s:e] += 1
    write(last_chrom, d, lens, prefix)
    plot_sample_counts(sample_counts, prefix)


def fmtr(x, p):
      v, suf = "%.2f" % (x / 1000000.), "M"
      # strip trailing 0's and "."
      while v[-1] in '.0' and '.' in v:
          v = v[:-1]
      return v + suf

def absfmtr(y, p):
    v = str(abs(y))
    while v[-1] in '.0' and '.' in v:
        v = v[:-1]
    return v

def plot_sample_counts(sample_counts, prefix):
    fig, ax = plt.subplots(1)
    ax.hist([np.array(sample_counts['1'].values()) / 2.0,
             np.array(sample_counts['2'].values()) / 2.0],
             20, label=['male', 'female'])
    ax.set_xlabel('Recombinations per meiosis')
    ax.set_ylabel('Count', rotation='vertical')
    plt.legend()
    plt.savefig('%s.recombinations-per-parent.png' % prefix)
    plt.close()


def write(chrom, d, lens, prefix):

    fig, ax = plt.subplots(1, figsize=(12, 3))

    for k, (sex, sex_label) in enumerate((('1', 'male'), ('2', 'female'), ('3', 'both'))):
        xs, ys = array.array('I'), array.array('I')
        if sex_label == 'both':
            # we output the total count for male+female
            # only print for the both and only plot for male, female separate
            arr = d['1'] + d['2']
        else:
            arr = d[sex]
        diff, = np.where(arr[:-1] != arr[1:])
        diff = list(diff + 1)
        for i, posn in enumerate(diff):
            if i == len(diff) - 1: break
            v = arr[posn]
            if v == 0: continue
            end = diff[i+1]
            vals = arr[posn:end]
            assert len(set(vals)) == 1, (vals, i, len(diff))
            if sex_label == "both":
                print("%s\t%d\t%d\t%d" % (chrom, posn, end, vals[0]))
            xs.extend((posn, end))
            ys.extend((vals[0], vals[0]))

        posn = diff[-1]
        if arr[posn] != 0:
            xs.extend((posn, posn + 1))
            ys.extend((arr[posn], arr[posn]))
            if sex_label == "both":
                print("%s\t%d\t%d\t%d" % (chrom, posn, posn + 1, arr[posn]))

        if sex_label == "both": break
        ys = np.asarray(ys, dtype=int)

        if k == 0:
            ys = -ys
        ax.plot(xs, ys, '-', label=sex_label)

    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(fmtr))
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(absfmtr))
    ax.legend(title="chromosome: " + chrom)
    ax.set_xlabel('Genomic position')
    ax.set_ylabel('Samples with crossover', rotation='vertical')

    plt.tight_layout()
    plt.savefig("%s%s.png" % (prefix.rstrip("."), chrom if prefix.endswith("/") else ("." + chrom)))
    plt.close(fig)

if __name__ == "__main__":
    import sys
    f = sys.argv[1]
    ped = sys.argv[2]
    prefix = sys.argv[3]
    main(f, ped, prefix)