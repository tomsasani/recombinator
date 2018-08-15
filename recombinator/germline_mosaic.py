from __future__ import print_function, absolute_import, division
import os
import sys
import time
from collections import OrderedDict

import pysam
from peddy import Ped
from cyvcf2 import VCF, Writer

import numpy as np
import scipy.stats as ss

def tranche99(filt, cutoff=99.6):
    """
    return True if the tranche is below 99.6
    VQSRTrancheINDEL90.00to99.00
    """
    if filt is None: return True
    if filt[:4] != "VQSR": return False
    try:
        return float(filt.split("to")[1]) < cutoff
    except:
        return False

def variant_prefilter(v, min_variant_qual):
    if len(v.REF) > 1: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if any([len(var) > 1 for var in v.ALT]): return False
    #if v.FILTER is not None and not tranche99(v.FILTER) : return False
    if v.QUAL < min_variant_qual: return False
    # if float(v.INFO.get('MQ')) < 30: return False
    return True

def validate_in_f2(v, sample_dict, f1, children, multiallelic=False):
    """
    if calling DNMs in the F1 generation, count the
    number of F2s that inherit that DNM

    v: cyvcf2 "Variant" object
    sample_dict: mapping of sample IDs to indices in the VCF FORMAT field 
    f1: F1, represented as peddy "sample" object
    children: list of peddy.Samples() for whom the F1 is a parent
    """

    # regardless of the DNM allele, below genotypes are constants
    gts = v.genotypes
    UNKNOWN = (-1, -1)
    HOM_REF = (0, 0)
    HET = (0, 1)
    HOM_ALT = (1, 1)

    if multiallelic:
        HET = (0, 2)
        HOM_ALT = (2, 2)

    quals = v.gt_quals

    f2_with_v = []
    for f2 in children:
        if f1.sample_id not in (f2.mom.sample_id, f2.dad.sample_id): continue
        f2_idx = sample_dict[f2.sample_id]
        mi, di = sample_dict[f2.mom.sample_id], sample_dict[f2.dad.sample_id]
        # redefine genotypes as tuples    
        m_gt, d_gt = (gts[mi][0], gts[mi][1]), (gts[di][0], gts[di][1])
        f2_gt = (gts[f2_idx][0], gts[f2_idx][1])
        if any([x == UNKNOWN for x in (m_gt, d_gt)]): continue
        # if we're on an autosome, the F2 should be HET for the DNM
        if v.CHROM != 'X' and f2_gt != HET: continue
        # if we're on the X, female F2s should be HET and male F2s should be HOM_ALT.
        elif v.CHROM == 'X':
            if f2.sex == 'male' and f2_gt != HOM_ALT: continue
            elif f2.sex == 'female' and f2_gt != HET: continue
        f2_with_v.append(f2)

    return f2_with_v

def get_denovo(v, sample_dict, p0, f1s, f2s,
        min_depth=5,
        min_qual=1,
        exclude=None,
        HET=1):
    """
    v: cyvcf2.Variant
    samples: dictionary of sample: index in vcf.
    f1s: list of peddy.Samples() that have parents, but not grandparents
    f2s: list of peddy.Samples() that have both parents and grandparents (we use these
    genotypes to corroborate de novo calls made in the F1s)
    min_depth: require all members of trio to have at least this depth.
    exclude: interval tree of regions to exclude.
    """

    if exclude is not None:
        if len(exclude[v.CHROM].search(v.start, v.end)) > 0:
            return None

    #if _use_cohort_filters and v.num_hom_alt > 1: return None

    ret = []
    gts = v.genotypes
   
    # when we filter out potential cohort DNMs, we want
    # to look in everyone but the grandparents
    all_children = f1s + f2s
    all_samples = p0 + f1s + f2s

    kids = f1s

    depths = None
    ref_depths = None
    # loop over alternate alleles
    for alt_index, k in enumerate(v.ALT):
        alt_depths = None
        # and then loop over the f1s. 
        for kid in kids:
            # Added check since sample '8477' apparently doesn't exist in our ceph data
            if kid.mom.sample_id == '8477' or kid.dad.sample_id == '8477':
                continue
            # basic sanity check for female f1s
            if v.CHROM == 'Y': continue

            UNKNOWN = (-1, -1)
            HOM_REF = (0, 0)
            HET = (0, 1)
            HOM_ALT = (1, 1)
            multiallelic = False
            if alt_index > 0: 
                multiallelic = True
                HET = (0, 2)
                HOM_ALT = (2, 2)

            quals = v.gt_quals
            quals[quals < 0] == 0

            if depths is None:
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0

            ki = sample_dict[kid.sample_id]

            # if we're looking in a third generation, kids (F1s) have children
            children = [k for k in f2s if kid.sample_id in (k.dad.sample_id, k.mom.sample_id)]

            if len(children) == 0: continue
            # need to check spouse for evidence of the mutation
            spouse = None
            for c in children:
                if kid == c.mom: spouse = c.dad
                elif kid == c.dad: spouse = c.mom
            si = sample_dict[spouse.sample_id]

            if spouse.mom is None or spouse.dad is None: continue

            sgmi, sgdi, kgmi, kgdi = (sample_dict[i] for i in (spouse.mom.sample_id,
                                                              spouse.dad.sample_id,
                                                              kid.mom.sample_id,
                                                              kid.dad.sample_id))

            progenitor_gts = [(gts[i][0], gts[i][1]) for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            progenitor_gqs = [quals[i] for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            progenitor_dps = [ref_depths[i] + alt_depths[i] for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            # minimal check on genotypes and qualities in the progenitors
            if not all([x == HOM_REF for x in progenitor_gts]): continue
            if not all([x >= min_qual for x in progenitor_gqs]): continue
            if not all([x >= min_depth for x in progenitor_dps]): continue

            # allele balance filter in progenitors
            grandparental_abs = []
            for i in (sgmi, sgdi, kgmi, kgdi):
                grandparental_abs.append((float(alt_depths[i]) / (alt_depths[i] + ref_depths[i])))

            k_ab = float(alt_depths[ki]) / (alt_depths[ki] + ref_depths[ki])
            s_ab = float(alt_depths[si]) / (alt_depths[si] + ref_depths[si])

            # loop over all samples outside of the family we're considering
            # and count carriers of the DNM allele
            possible_carriers = 0 # genotype only
            likely_carriers = 0 # must meet QUAL and AB threshold
            
            for sample in all_samples:
                if sample.family_id == kid.family_id: continue
                sample_idx = sample_dict[sample.sample_id]
                genotype = (gts[sample_idx][0], gts[sample_idx][1])
                if genotype != HET: continue
                if v.CHROM == 'X':
                    if sample.sex == 'male' and genotype != HOM_ALT: continue
                possible_carriers += 1
                sample_ref = ref_depths[sample_idx]
                sample_alt = alt_depths[sample_idx]
                sample_total = float(sample_ref + sample_alt)
                sample_qual = quals[sample_idx]
                # "likely" carriers meet stricter threshold
                if sample_qual < 20: continue
                if sample_total < 10: continue
                if v.CHROM == 'X':
                    if sample.sex == 'male' and sample_alt / sample_total < 0.75: continue
                else:
                    if not 0.25 <= (sample_alt / sample_total) <= 0.75: continue
                likely_carriers += 1
            possible_carriers -= likely_carriers

            # check for presence of the DNM in a third generation
            f2_with_v = []
            f2_with_v = validate_in_f2(v, sample_dict, kid, children, multiallelic=multiallelic)
            if len(f2_with_v) < 1: continue
            f2_with_v_ab = []
            f2_with_v_gq = []
            f2_with_v_dp = []
            # track allele balance and qualities in F2 DNM events
            for f2 in f2_with_v:
                f2_idx = sample_dict[f2.sample_id]
                total_depth = alt_depths[f2_idx] + ref_depths[f2_idx]
                f2_with_v_dp.append(total_depth)
                ab = alt_depths[f2_idx] / float(total_depth)
                gq = quals[f2_idx]
                f2_with_v_ab.append(ab)
                f2_with_v_gq.append(gq)
            ret.extend(variant_info(v, kid, sample_dict, f2_with_v=f2_with_v, 
                                                        f2_with_v_ab=f2_with_v_ab,
                                                        f2_with_v_gq=f2_with_v_gq,
                                                        f2_with_v_dp=f2_with_v_dp,
                                                        alt_i=alt_index,
                                                        grandparental_abs=grandparental_abs,
                                                        k_ab=k_ab, s_ab=s_ab,
                                                        possible_carriers=possible_carriers,
                                                        likely_carriers=likely_carriers))
            
    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    #print (len(ret))
    if len(ret) > 0: 
        return ret

def variant_info(v, kid, sample_dict, f2_with_v=[], f2_with_v_ab=[], f2_with_v_gq=[], f2_with_v_dp=[],
                    interval_type=None, alt_i=None,
                    grandparental_abs=[], k_ab=0., s_ab=0., possible_carriers=0, likely_carriers=0):

    ki, mi, di = sample_dict[kid.sample_id], sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

    f2_with_v_idx = [sample_dict[f2.sample_id] for f2 in f2_with_v]

    ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals

    k_ref, k_alt = ref_depths[ki], alt_depths[ki]

    ac = v.INFO.get('AC')
    if type(ac) is tuple:
        ac = ac[alt_i]

    yield OrderedDict((
        ("chrom", v.CHROM),
        ("start", v.start),
        ("end", v.end),
        ("sample_id", kid.sample_id),
        ("family_id", kid.family_id),
        ("sample_sex", kid.sex),
        ("ref", v.REF),
        ("alt", v.ALT[alt_i]),
        ("alt_i", str(alt_i + 1) + '/' + str(len(v.ALT))),
        ("num_hom_alt", v.num_hom_alt),
        ("ac", ac),
        ("an", v.INFO.get("AN")),
        ("possible_carriers", possible_carriers),
        ("likely_carriers", likely_carriers),
        ("mq_vcf", v.INFO.get("MQ")),
        ("mapq_ranksum_vcf", v.INFO.get('MQRankSum')),
        ("base_q_ranksum_vcf", v.INFO.get("BaseQRankSum")),
        ("clipping_ranksum_vcf", v.INFO.get("ClippingRankSum")),
        ("excess_het_vcf", v.INFO.get("ExcessHet")),
        ("strand_bias_vcf", v.INFO.get("FS")),
        ("inbreeding_coef_vcf", v.INFO.get("InbreedingCoeff")),
        ("qual_by_depth", v.INFO.get("QD")),
        ("raw_mq_vcf", v.INFO.get("RAW_MQ")),
        ("read_pos_ranksum_vcf", v.INFO.get("ReadPosRankSum")),
        ("symmetric_odd_ratio_sb_vcf", v.INFO.get("SOR")),
        ("vqslod_vcf", v.INFO.get("VQSLOD")),
        ("neg_train_site", v.INFO.get("NEGATIVE_TRAIN_SITE")),
        ("pos_train_site", v.INFO.get("POSITIVE_TRAIN_SITE")),
        ("filter", v.FILTER), #or "PASS"),
        ("paternal_id", kid.paternal_id),
        ("maternal_id", kid.maternal_id),
        ("kid_ref_depth", k_ref),
        ("kid_alt_depth", k_alt),
        ("kid_total_depth", k_ref + k_alt),
        ("kid_allele_balance", k_alt / float(k_alt + k_ref)),
        ("mom_ref_depth", ref_depths[mi]),
        ("mom_alt_depth", alt_depths[mi]),
        ("mom_total_depth", ref_depths[mi] + alt_depths[mi]),
        ("dad_ref_depth", ref_depths[di]),
        ("dad_alt_depth", alt_depths[di]),
        ("dad_total_depth", ref_depths[di] + alt_depths[di]),
        ("kid_qual", quals[ki]),
        ("mom_qual", quals[mi]),
        ("dad_qual", quals[di]),
        ("depth_mean", "%.1f" % np.mean(ref_depths + alt_depths)),
        ("depth_std", "%.1f" % np.std(ref_depths + alt_depths)),
        ("qual_mean", "%.1f" % np.mean(quals)),
        ("call_rate", "%.3f" % v.call_rate),
        ("f2_with_variant", ','.join([str(x.sample_id) for x in f2_with_v])),
        ("f2_with_variant_ab", ','.join([str(x) for x in f2_with_v_ab])),
        ("f2_with_variant_dp", ','.join([str(x) for x in f2_with_v_dp])),
        ("f2_with_variant_gq", ','.join([str(x) for x in f2_with_v_gq])),
        ("grandparental_abs", ','.join([str(x) for x in grandparental_abs])),
        ("kid_ab", k_ab),
        ("spouse_ab", s_ab)
        ))

def denovo(v, sample_dict, p0, f1s, f2s,
        min_depth=10, 
        min_qual=1,
        exclude=None,
        HET=1):
    if not variant_prefilter(v, 1): return None
    return get_denovo(v, sample_dict, p0, f1s, f2s,
            min_depth=min_depth,
            min_qual=min_qual,
            exclude=exclude,
            HET=HET)


def write_denovo(d, fh, _line={}):
    if _line.get(fh.name, 0) == 0:
        _line[fh.name] = 1
        print("#" + "\t".join(d[0].keys()), file=fh)
    for x in d:
        print("\t".join(map(str, x.values())), file=fh)

def run(args):

    vcf = VCF(args.vcf)
    ped = Ped(args.ped)
    psamples = set([x.sample_id for x in ped.samples()])
    samples = [x for x in vcf.samples if x in psamples]
    if args.exclude:
        from .recombinator import read_exclude
        exclude = read_exclude(args.exclude, args.chrom)
    else:
        exclude = None

    vcf = VCF(args.vcf, gts012=True, samples=samples)
    wtr = Writer("-", vcf)

    sample_dict = {v: i for i, v in enumerate(vcf.samples)}
    kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    p0 = [k for k in ped.samples() if k.mom is None and k.dad is None]
    f1s = [k for k in kids if k.mom.mom is None and k.dad.dad is None and k.mom.dad is None and k.dad.mom is None]
    f2s = [k for k in kids if k not in f1s]

    _use_cohort_filters = False if os.environ.get("DN_NO_COHORT") else True

    n_dn = 0
    for i, v in enumerate(vcf(str(args.chrom)) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)

        d = denovo(v, sample_dict, p0, f1s, f2s, min_depth=args.min_depth,
                min_qual=args.min_qual,
                exclude=exclude)
        if d is not None:
            write_denovo(d, sys.stdout)
            n_dn += 1

    sys.stdout.flush()
    print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)


def main(argv):
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--chrom", help="extract only this chromosome", default=None)
    p.add_argument("-exclude", help="BED file of regions to exclude (e.g. LCRs)")
    p.add_argument("-min_depth", help="min depth in kid and both parents [10]", type=int, default=10)
    p.add_argument("-min_qual", help="minimum GQ in kid [0]", type=int, default=0)
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
