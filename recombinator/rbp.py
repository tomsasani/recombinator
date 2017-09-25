from __future__ import absolute_import, print_function
import sys
import pysam
import toolshed as ts
from cyvcf2 import VCF
import scipy.stats as ss
from .var_info import get_position


def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--fragment-length", type=int, default=350, help="library fragment length")
    p.add_argument("--ped-file", required=True, help="optional ped file if parents aren't specified in denovo file")
    p.add_argument("bed", help="bed file of denovos with a column header for 'chrom,start,end,sample_id")
    p.add_argument("vcf")
    args = p.parse_args(argv)
    run(args)

def phased(v, d, kid_id, ev_min=3, prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):

    evidence = 0

    v_start = int(v.POS) - 1 # The putative HET variant in mom/dad we're interested in phasing with...
    b_start = int(d['start']) # ...the DNM we've already identified as HET in the kid.
    dist = v_start - b_start

    bam = pysam.AlignmentFile(prefix + kid_id + ".bam", "rb")

    for pileupcolumn in bam.pileup(v.CHROM, min(v_start, b_start),
            max(v_start, b_start)):
        if pileupcolumn.pos != v_start:
            continue
        # For every read that spans the putative HET to phase,
        # check that the read has support for both the putative HET
        # AND the known DNM.
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            query = pileupread.alignment.query_sequence
            # If variant to phase is before the DNM.
            if dist < 0:
                try:
                    if query[pileupread.query_position].upper() == v.ALT[0].upper() and query[pileupread.query_position + abs(dist)] == d['alt']:
                        evidence += 1
                except IndexError:
                    continue
            else:
                try:
                    if query[pileupread.query_position].upper() == v.ALT[0].upper() and query[pileupread.query_position - abs(dist)] == d['alt']:
                        evidence += 1
                except IndexError:
                    continue

    if evidence >= ev_min:
        return True

def filter_sites(v, i, min_depth=10, min_ab_p=0.05):
    ref_depths, alt_depths = v.gt_ref_depths, v.gt_alt_depths
    pab = ss.binom_test([ref_depths[i], alt_depths[i]], p=0.5, alternative='two-sided')
    if ref_depths[i] + alt_depths[i] < min_depth:
        return None
    if pab < min_ab_p:
        return None 
    return True

def run(args):

    HET, HOM_REF = 1, 0

    vcf = VCF(args.vcf, gts012=True)
    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    sample_to_dad = {toks[1]: sample_lookup.get(toks[2]) for toks in ts.reader(args.ped_file, header=False)}
    sample_to_mom = {toks[1]: sample_lookup.get(toks[3]) for toks in ts.reader(args.ped_file, header=False)}

    for i, d in enumerate(ts.reader(args.bed, header="ordered")):

        idad = sample_to_dad[d['sample_id']]
        imom = sample_to_mom[d['sample_id']]
        ikid = sample_lookup[d['sample_id']]
        for v in get_position(vcf, d, extra=args.fragment_length):
            if len(v.ALT) > 1:
                continue
            gt_types = v.gt_types
            if gt_types[ikid] != HET:
                continue
            if gt_types[idad] == HET and gt_types[imom] == HOM_REF:
                #if filter_sites(v, idad) is not None and phased(v, d, d['sample_id']):
                if phased(v, d, d['sample_id']):
                    #d['dad-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))
                    d['phased_parent'] = d['paternal_id']
            elif gt_types[imom] == HET and gt_types[idad] == HOM_REF:
                #if filter_sites(v, imom) is not None and phased(v, d, d['sample_id']):
                if phased(v, d, d['sample_id']):
                    #d['mom-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))
                    d['phased_parent'] = d['maternal_id']
            else:
                d['phased_parent'] = '' 
        if i == 0:
            print("\t".join(d.keys()))
        #d['phased_parent'] = '\t'.join(d['phased_parent'])
        print("\t".join(d.values()))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    main(sys.argv[1:])
