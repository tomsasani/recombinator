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
    if len(v.REF) > 15: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if len(v.ALT[0]) > 15: return False
    if v.FILTER is not None and not tranche99(v.FILTER) : return False
    if v.QUAL < min_variant_qual: return False
    if float(v.INFO.get('MQ')) < 30: return False
    if type(v.INFO.get('AC')) is int:
        if v.INFO.get('AC') > 300: return False
    else:
        if min(v.INFO.get('AC')) > 300: return False
    return True

def validate_in_f2(v, sample_dict, gts, kid, f2s, multiallelic=False):
    """
    at this point, the variant is a solid putative DNM. now need
    to validate that it's in at least one of the f2s. all of the
    previous filters were on a population basis, so for now we'll only
    check if one of the kids is a HET
    """

    UNKNOWN = (-1, -1)
    HOM_REF = (0, 0)
    HET = (0, 1)
    HOM_ALT = (1, 1)

    if multiallelic:
        HET = (0, 2)
        HOM_ALT = (2, 2)

    quals = v.gt_quals
    f2_with_v = []
    for f2 in f2s:
        f2_idx = sample_dict[f2.sample_id]
        mi, di = sample_dict[f2.mom.sample_id], sample_dict[f2.dad.sample_id]

        m_gt, d_gt = (gts[mi][0], gts[mi][1]), (gts[di][0], gts[di][1])
        f2_gt = (gts[f2_idx][0], gts[f2_idx][1])
        if any([x == UNKNOWN for x in (m_gt, d_gt)]): continue
        if kid.sample_id not in [f2.mom.sample_id, f2.dad.sample_id]: continue

        # Make sure that the F1's partner doesn't also have the same DNM
        # Actually, this is possible? But extreeeeemely unlikely. This check only
        # works if we don't look at the X.
        if v.CHROM != 'X':
            # only allow for one DNM across both parents
            if multiallelic and sum([sum(x) for x in (m_gt, d_gt)]) > 2: continue
            elif not multiallelic and sum([sum(x) for x in (m_gt, d_gt)]) > 1: continue
            if f2_gt != HET: continue
        # Special logic for inherited DNM on male X.
        elif v.CHROM == 'X':
            # Female F1 with DNM passes it on to male F2, must be
            # homozygous in the F2.
            if kid.sex == 'female' and f2.sex == 'male':
                if f2_gt != HOM_ALT: continue
            # Male F1 with DNM on X simply cannot pass it onto male F2.
            elif kid.sex == 'male' and f2.sex == 'male': continue
            else:
                if f2_gt != HET: continue

        f2_with_v.append(f2)

    return f2_with_v

def get_ab(samp_id, ref_depths, alt_depths, sample_dict):
    idx = sample_dict[samp_id]
    ab = float(alt_depths[idx]) / (alt_depths[idx] + ref_depths[idx])
    return ab


def get_denovo(v, sample_dict, p0, f1s, f2s,
        max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_qual=1,
        exclude=None,
        _use_cohort_filters=True,
        simons_prefix=False,
        check_bams=False,
        HET=1):
    """
    v: cyvcf2.Variant
    samples: dictionary of sample: index in vcf.
    f1s: list of peddy.Samples() that have parents, but not grandparents
    f2s: list of peddy.Samples() that have both parents and grandparents (we use these
    genotypes to corroborate de novo calls made in the F1s)
    max_alts_in_parents: max number of alternate reads that can appear in the parents.
    min_depth: require all members of trio to have at least this depth.
    min_allele_balance_p: if the p-value for ref, alt counts is less than this, exclude
    min_depth_percentile: require all members of a trio to be in this percentile of depth.
    exclude: interval tree of regions to exclude.
    _use_cohort_filters: if set to false then only use info within a family to
    decide if a variant is denovo. don't use things like AF in the cohort.
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
            if simons_prefix and 'SSC' not in kid.sample_id: continue
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

            parental_gts = [(gts[i][0], gts[i][1]) for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            parental_gqs = [quals[i] for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            parental_dps = [ref_depths[i] + alt_depths[i] for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]
            if not all([x == HOM_REF for x in parental_gts]): continue
            if not all([x > min_qual for x in parental_gqs]): continue
            if not all([x > min_depth for x in parental_dps]): continue

            # allele balance filter in progenitors
            parental_abs = [(float(alt_depths[i]) / (alt_depths[i] + ref_depths[i])) for i in (ki, si, sgmi, sgdi, kgmi, kgdi)]

            family_idxs = [sample_dict[k.sample_id] for k in all_samples if k.family_id == kid.family_id]
            gt_to_filter = HET
            if v.CHROM == 'X' and kid.sex == 'male': gt_to_filter = HOM_ALT
            possible_carriers = 0
            likely_carriers = 0

            for sample in all_samples:
                if sample.family_id == kid.family_id: continue
                sample_idx = sample_dict[sample.sample_id]
                genotype = (gts[sample_idx][0], gts[sample_idx][1])
                if genotype != gt_to_filter: continue
                possible_carriers += 1
                sample_ref = ref_depths[sample_idx]
                sample_alt = alt_depths[sample_idx]
                sample_total = float(sample_alt + sample_ref)
                sample_qual = quals[sample_idx]

                if sample_qual < 10: continue
                if sample_total < 10: continue
                if gt_to_filter == HOM_ALT:
                    if sample_alt / sample_total < 0.75: continue
                elif gt_to_filter == HET:
                    if not 0.25 <= (sample_alt / sample_total) <= 0.75: continue
                likely_carriers += 1
            possible_carriers -= likely_carriers
            
            # check for inheritance of the DNM in a third generation
            f2_with_v = []
            f2_with_v = validate_in_f2(v, sample_dict, gts, kid, children, multiallelic=multiallelic)
            if len(f2_with_v) < 2: continue
            f2_with_v_ab = []
            f2_with_v_gq = []
            f2_with_v_dp = []
            # impose lenient filters on individuals that inherit the DNM
            # track allele balance in transmission events to look for mosaics
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
                                                        possible_carriers=possible_carriers,
                                                        likely_carriers=likely_carriers,
                                                        k_gm_ab=parental_abs[4], k_gd_ab=parental_abs[5], 
                                                        s_gm_ab=parental_abs[2], s_gd_ab=parental_abs[3], 
                                                        k_ab=parental_abs[0], s_ab=parental_abs[1]))
            
    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    #print (len(ret))
    if len(ret) > 0: 
        return ret

def variant_info(v, kid, sample_dict, f2_with_v=[], f2_with_v_ab=[], f2_with_v_gq=[], f2_with_v_dp=[],
                    interval_type=None, alt_i=None, p_dnm_in_f2s=0.,
                    sb_p=0., possible_carriers=0., likely_carriers=0.,
                    k_gm_ab=0., k_gd_ab=0., s_gm_ab=0., s_gd_ab=0., k_ab=0., s_ab=0.):

    ki, mi, di = sample_dict[kid.sample_id], sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

    f2_with_v_idx = [sample_dict[f2.sample_id] for f2 in f2_with_v]

    ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals

    k_ref, k_alt = ref_depths[ki], alt_depths[ki]

    alt_sum = sum([alt_depths[i] for i in range(len(alt_depths)) if i not in f2_with_v_idx]) - k_alt
    ref_sum = sum([ref_depths[i] for i in range(len(ref_depths)) if i not in f2_with_v_idx]) - k_ref

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
        ("possible_carriers", possible_carriers),
        ("likely_carriers", likely_carriers),
        ("num_hom_alt", v.num_hom_alt),
        ("allele_number", v.INFO.get("AN")),
        ("mq_vcf", v.INFO.get("MQ")),
        ("mapq_ranksum_vcf", v.INFO.get('MQRankSum')),
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
        ("cohort_alt_depth", alt_sum),
        ("f2_with_variant", ','.join([str(x.sample_id) for x in f2_with_v])),
        ("f2_with_variant_ab", ','.join([str(x) for x in f2_with_v_ab])),
        ("f2_with_variant_dp", ','.join([str(x) for x in f2_with_v_dp])),
        ("f2_with_variant_gq", ','.join([str(x) for x in f2_with_v_gq])),
        ("kid_grandma_ab", k_gm_ab),
        ("kid_grandpa_ab", k_gd_ab),
        ("spouse_grandma_ab", s_gm_ab),
        ("spouse_grandpa_ab", s_gd_ab),
        ("kid_ab", k_ab),
        ("spouse_ab", s_ab)
        ))

def denovo(v, sample_dict, p0, f1s, f2s,
        min_depth=10, 
        max_mean_depth=400,
        min_qual=1,
        _use_cohort_filters=True,
        exclude=None,
        simons_prefix=False,
        check_bams=False,
        HET=1):
    if not variant_prefilter(v, 1): return None
    return get_denovo(v, sample_dict, p0, f1s, f2s,
            min_depth=min_depth,
            max_mean_depth=max_mean_depth,
            min_qual=min_qual,
            exclude=exclude,
            _use_cohort_filters=_use_cohort_filters,
            simons_prefix=simons_prefix,
            check_bams=check_bams,
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

    simons_prefix, check_bams = False, False
    if args.simons_prefix:
        simons_prefix = True
    if args.check_bams:
        check_bams = True
    n_dn = 0
    for i, v in enumerate(vcf(str(args.chrom)) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)

        d = denovo(v, sample_dict, p0, f1s, f2s, min_depth=args.min_depth,
                min_qual=args.min_qual,
                exclude=exclude, _use_cohort_filters=_use_cohort_filters,
                simons_prefix=simons_prefix, check_bams=check_bams)
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
    p.add_argument("-min_qual", help="minimum GQ in kid [1]", type=int, default=1)
    p.add_argument("-simons_prefix", help="require all sample names to begin with the 'SSC' prefix", action="store_true")
    p.add_argument("-check_bams", help="check parental BAMs for evidence of DNM", action="store_true")
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
