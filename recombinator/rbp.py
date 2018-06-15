from __future__ import absolute_import, print_function
import sys
import pysam
import toolshed as ts
from cyvcf2 import VCF
from peddy import Ped
from collections import defaultdict
import scipy.stats as ss

def main(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("-fragment-length", type=int, default=500, help="library fragment length")
    p.add_argument("-chrom")
    p.add_argument("-prefix", help='path to directory containing BAM files, formatted as <sample_name.bam>', default="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/")
    p.add_argument("bed", help="bed file of denovos with a column header for 'chrom,start,end,sample_id")
    p.add_argument("vcf", help="VCF file containing (at least) all variants for the trio samples")
    args = p.parse_args(argv)
    run(args)

def get_position(vcf, d, ref=None, alt=None, extra=0):
    loc = "%s:%d-%d" % (d['chrom'], int(d['start']) + 1 - extra,
                        int(d['end']) + extra)
    for v in vcf(loc):
        if ref and v.REF != ref: continue
        if alt and alt not in v.ALT: continue
        yield v

def calculate_evidence(v, d, inf_base='A', dnm_base='T'):

    parent_evidence, other_parent_evidence = 0, 0

    if inf_base.upper() == v.ALT[0].upper() and dnm_base == d['alt']:
        parent_evidence += 1
    # if the read supports the reference at both sites, we know that other reads with
    # the DNM are phased with the alternate allele.
    elif inf_base.upper() == v.REF.upper() and dnm_base == d['ref']:
        parent_evidence += 1
    # NOTE: evidence for other parent's chromosome, if
    # the read supports the DNM (d[alt]) but has reference support at the
    # putative HET variant
    elif inf_base.upper() == v.REF.upper() and dnm_base == d['alt']:
        other_parent_evidence += 1

    return parent_evidence, other_parent_evidence

def high_quality_site(v, i, sex='male', min_depth=10, min_ab=0.25):
    min_ab = 0.75 if v.CHROM == 'X' and sex == 'male' else min_ab
    ref_depths, alt_depths = v.gt_ref_depths, v.gt_alt_depths
    if alt_depths[i] / float(ref_depths[i] + alt_depths[i]) < min_ab: 
        return False
    if ref_depths[i] + alt_depths[i] < min_depth:
        return False
    return True

def phased(v, d, kid_id, prefix="./"):

    read_names = defaultdict() 

    parent_evidence, other_parent_evidence = 0, 0

    # The informative variant in mom/dad we're interested in phasing with...
    inf_ref_pos = v.start
    # ...the DNM we've already identified as HET in the kid.
    dnm_ref_pos = int(d['start']) 
    dist = inf_ref_pos - dnm_ref_pos

    bam = pysam.AlignmentFile(prefix + kid_id + ".bam", "rb")

    # loop over pileups in the BAM until we hit the position of the putative HET to phase with
    for pileupcolumn in bam.pileup(v.CHROM, min(inf_ref_pos, dnm_ref_pos),
            max(inf_ref_pos, dnm_ref_pos)):
        # do everything w/r/t the informative site
        if pileupcolumn.pos != inf_ref_pos: continue
        # for every read that spans the informative site,
        # check that the read has support for both the putative HET
        # AND the known DNM.
        for pileupread in pileupcolumn.pileups:
            # ignore reads with indels w/r/t the reference
            if pileupread.is_del or pileupread.is_refskip: continue
            if pileupread.alignment.mapping_quality < 20: continue
            # get the alignment's sequence
            l_query_seq = pileupread.alignment.query_sequence
            l_start, l_end = pileupread.alignment.reference_start, pileupread.alignment.reference_end

            inf_query_pos = pileupread.query_position
            inf_query_base = l_query_seq[inf_query_pos]

            # if only one of the two variant sites are in the same read
            if dnm_ref_pos not in range(l_start, l_end):
                try: mate = bam.mate(pileupread.alignment)
                # if we can't find the mate, we can't phase this variant
                except ValueError: continue
                r_query_seq = mate.query_sequence
                r_start, r_end = mate.reference_start, mate.reference_end
                # if the DNM isn't in the mate sequence, either, skip
                if dnm_ref_pos not in range(r_start, r_end): continue
                dnm_query_pos = dnm_ref_pos - r_start
                try:
                    dnm_query_base = r_query_seq[dnm_query_pos]
                except IndexError: continue

                p, o = calculate_evidence(v, d, inf_base=inf_query_base, dnm_base=dnm_query_base)
                parent_evidence += p
                other_parent_evidence += o
            else:
                dnm_query_pos = inf_query_pos - dist
                try:
                    dnm_query_base = l_query_seq[dnm_query_pos]
                except IndexError: continue
                p, o = calculate_evidence(v, d, inf_base=inf_query_base, dnm_base=dnm_query_base)
                parent_evidence += p
                other_parent_evidence += o

    return parent_evidence, other_parent_evidence

def run(args):

    HOM_ALT, HET, HOM_REF = 3, 1, 0

    vcf = VCF(args.vcf)
    sample_lookup = {s: i for i, s in enumerate(vcf.samples)}

    for i, d in enumerate(ts.reader(args.bed, header="ordered")):

        if args.chrom and d['chrom'] != args.chrom: continue

        ikid = sample_lookup[d['sample_id']]
        idad = sample_lookup[d['paternal_id']]
        imom = sample_lookup[d['maternal_id']]

        # loop over all variants in the VCF +/- 350 bp from the DNM
        phased_parent, phased_variant = '', ''
        evidence = '' 
        for v in get_position(vcf, d, extra=args.fragment_length):
            # ignore more complex variants for now
            if len(v.REF) > 1: continue
            if len(v.ALT[0]) > 1: continue
            gt_types = v.gt_types
            # special logic for male chrX (hemizygous) variants
            is_male = d['sample_sex'] == 'male'
            if v.CHROM == 'X' and is_male and gt_types[ikid] != HOM_ALT: continue
            else:
                if gt_types[ikid] != HET: continue

            if gt_types[idad] in (HET, HOM_ALT) and gt_types[imom] == HOM_REF:
                if not high_quality_site(v, idad): continue
                dad_evidence, mom_evidence = phased(v, d, d['sample_id'], prefix=args.prefix)
                if mom_evidence + dad_evidence < 3: continue
                if abs(mom_evidence - dad_evidence) != max([mom_evidence, dad_evidence]): continue
                if dad_evidence > mom_evidence: 
                    phased_parent = d['paternal_id']
                # assume that if no evidence on paternal ID (and it's a HET), must
                # be from mom.
                # TODO: check that this is true!!!
                elif dad_evidence < mom_evidence:
                    phased_parent = d['maternal_id']
                else: continue
                evidence = str(dad_evidence) + ',' + str(mom_evidence)
                phased_variant = str(v.start)

            elif gt_types[imom] in (HET, HOM_ALT) and gt_types[idad] == HOM_REF:
                if not high_quality_site(v, imom, sex='female'): continue
                mom_evidence, dad_evidence = phased(v, d, d['sample_id'], prefix=args.prefix)
                if mom_evidence + dad_evidence < 3: continue
                if abs(mom_evidence - dad_evidence) != max([mom_evidence, dad_evidence]): continue
                if mom_evidence > dad_evidence:
                    phased_parent = d['maternal_id']
                elif mom_evidence < dad_evidence:
                    phased_parent = d['paternal_id']
                else: continue
                evidence = str(dad_evidence) + ',' + str(mom_evidence)
                phased_variant = str(v.start)

        d['phased_parent'] = phased_parent
        d['phased_variant'] = phased_variant
        d['evidence_count_p,m'] = evidence
        if i == 0:
            print("\t".join(d.keys()))
        print("\t".join(d.values()))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    main(sys.argv[1:])
