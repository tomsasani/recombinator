from __future__ import absolute_import, print_function
import sys
import pysam
import toolshed as ts
from cyvcf2 import VCF
from collections import defaultdict
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

def phased(v, d, kid_id, het=True, prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):

    read_names = defaultdict() 

    parent_evidence, other_parent_evidence = 0, 0

    rcdict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}

    het_reference_pos = int(v.POS) - 1 # The putative HET variant in mom/dad we're interested in phasing with...
    dnm_reference_pos = int(d['start']) # ...the DNM we've already identified as HET in the kid.
    dist = het_reference_pos - dnm_reference_pos

    #search_mate = False
    #if abs(dist) > 151: search_mate = True

    bam = pysam.AlignmentFile(prefix + kid_id + ".bam", "rb")

    for pileupcolumn in bam.pileup(v.CHROM, min(het_reference_pos, dnm_reference_pos),
            max(het_reference_pos, dnm_reference_pos)):
        if pileupcolumn.pos != het_reference_pos:
            continue
        # For every read that spans the putative HET to phase,
        # check that the read has support for both the putative HET
        # AND the known DNM.
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            l_query_seq = pileupread.alignment.query_sequence
            l_start, l_end = pileupread.alignment.reference_start, pileupread.alignment.reference_end
            try:
                mate = bam.mate(pileupread.alignment)
            except ValueError: continue
            r_query_seq = mate.query_sequence
            r_start, r_end = mate.reference_start, mate.reference_end

            l_range = list(range(l_start, l_end))
            r_range = list(range(r_start, r_end))
            insert = max(l_start, r_start) - min(l_end, r_end)

            if dnm_reference_pos not in l_range and dnm_reference_pos not in r_range: continue
           
            het_query_pos = pileupread.query_position
            # position of known DNM
            dnm_query_pos = het_query_pos - dist
            # query base at putative HET variant
            query_het_base = l_query_seq[het_query_pos]
            if dnm_reference_pos in r_range:
                # query position in the mate of the DNM
                try:
                    dnm_query_pos = dnm_reference_pos - r_start
                    query_dnm_base = r_query_seq[dnm_query_pos]
                except IndexError:
                    continue
                    #print (r_range, file=sys.stderr)
                    #print ("dnm reference pos", dnm_reference_pos, file=sys.stderr)
            elif dnm_reference_pos in l_range:
                # query base at known DNM
                try:
                    query_dnm_base = l_query_seq[dnm_query_pos]
                except IndexError:
                    continue
            if query_het_base.upper() == v.ALT[0].upper() and query_dnm_base == d['alt']:
                parent_evidence += 1
                # NOTE: evidence for other parent's chromosome, if
                # the read supports the DNM (d[alt]) but has reference support at the
                # putative HET variant
            if het and query_het_base.upper() == v.REF[0].upper() and query_dnm_base == d['alt']:
                    other_parent_evidence += 1

    return parent_evidence, other_parent_evidence

def filter_sites(v, i, min_depth=10, min_ab=0.1):
    ref_depths, alt_depths = v.gt_ref_depths, v.gt_alt_depths
    if alt_depths[i] / float(ref_depths[i] + alt_depths[i]) < min_ab: 
        return None
    if ref_depths[i] + alt_depths[i] < min_depth:
        return None
    return True

def run(args):

    HOM_ALT, HET, HOM_REF = 3, 1, 0

    vcf = VCF(args.vcf)
    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    sample_to_dad = {toks[1]: sample_lookup.get(toks[2]) for toks in ts.reader(args.ped_file, header=False)}
    sample_to_mom = {toks[1]: sample_lookup.get(toks[3]) for toks in ts.reader(args.ped_file, header=False)}

    for i, d in enumerate(ts.reader(args.bed, header="ordered")):

        idad = sample_to_dad[d['sample_id']]
        imom = sample_to_mom[d['sample_id']]
        ikid = sample_lookup[d['sample_id']]
        # loop over all variants in the VCF +/- 350 bp from the 
        # known DNM
        for v in get_position(vcf, d, extra=args.fragment_length):
            if len(v.REF) > 1: continue
            if len(v.ALT[0]) > 1: continue
            gt_types = v.gt_types
            if gt_types[ikid] != HET:
                continue
            if gt_types[idad] in (HET, HOM_ALT) and gt_types[imom] == HOM_REF:
                #if filter_sites(v, idad) is not None and phased(v, d, d['sample_id']):
                dad_evidence, mom_evidence = phased(v, d, d['sample_id'])
                    #d['dad-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))
                if dad_evidence > mom_evidence:
                    d['phased_parent'] = d['paternal_id']
                # assume that if no evidence on paternal ID (and it's a HET), must
                # be from mom.
                # TODO: check that this is true!!!
                elif dad_evidence < mom_evidence:
                    d['phased_parent'] = d['maternal_id']
                else: continue
                d['phased_variant'] = str(v.start)
            elif gt_types[imom] in (HET, HOM_ALT) and gt_types[idad] == HOM_REF:
                #if filter_sites(v, imom) is not None and phased(v, d, d['sample_id']):
                if mom_evidence > dad_evidence:
                    #d['mom-sites'].append("%s:%d-%d" % (v.CHROM, v.start+1, v.end))
                    d['phased_parent'] = d['maternal_id']
                elif mom_evidence < dad_evidence:
                    d['phased_parent'] = d['paternal_id']
                else: continue
                d['phased_variant'] = str(v.start)
            else:
                d['phased_parent'] = '' 
                d['phased_variant'] = ''
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
