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
    if v.FILTER is not None and not tranche99(v.FILTER) : return False

    if v.QUAL < min_variant_qual: return False
    return True


def get_ab(rd, ad, idx):
    return ad[idx] / float(rd[idx] + ad[idx])

def get_denovo(v, sample_dict, p0, f1s, f2s,
        max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_qual=1,
        multigen=False,
        exclude=None,
        _use_cohort_filters=True,
        simons_prefix=False,
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

    HOM_REF, HET, HOM_ALT = 0, 1, 2

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
            ki = sample_dict[kid.sample_id]
            mi, di = sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]
            ref_depths, alt_depths = v.gt_ref_depths, v.gt_alt_depths

            quals = v.gt_quals
            quals[quals < 0] == 0
            # basic sanity check for female f1s
            if v.CHROM == 'Y': continue

            # define genotypes
            UNKNOWN = (-1, -1)
            HOM_REF = (0, 0)
            HET = (0, 1)
            HOM_ALT = (1, 1) 
            
            multiallelic = False
            if alt_index > 0: 
                multiallelic = True
                HET = (0, 2)
                HOM_ALT = (2, 2)
            
            # F1 must have children
            children = [k for k in f2s if kid.sample_id in (k.dad.sample_id, k.mom.sample_id)]
            if len(children) == 0: continue

            # get genotypes of proband and their parents
            k_gt, m_gt, d_gt = (gts[ki][0], gts[ki][1]),\
                               (gts[mi][0], gts[mi][1]),\
                               (gts[di][0], gts[di][1])
           
            # find the proband's spouse
            spouse = None
            for c in children:
                if kid == c.mom: spouse = c.dad
                elif kid == c.dad: spouse = c.mom
            si = sample_dict[spouse.sample_id]
           
            # make sure the proband's spouse has parents in the cohort
            if spouse.mom is None or spouse.dad is None: continue

            smi, sdi = sample_dict[spouse.mom.sample_id], sample_dict[spouse.dad.sample_id]
           
            # get genotypes of proband's spouse and their parents
            s_gt, sm_gt, sd_gt = (gts[si][0], gts[si][1]),\
                               (gts[smi][0], gts[smi][1]),\
                               (gts[sdi][0], gts[sdi][1])

            # skip sites where any of the queried individuals have UNKNOWN genotypes
            if any([gt == UNKNOWN for gt in (k_gt, m_gt, d_gt,
                                             s_gt, sm_gt, sd_gt)]): continue

            # skip sites where the F1 individuals are not HOM_REF
            if not (k_gt == HOM_REF and s_gt == HOM_REF): continue

            quals = v.gt_quals
            quals[quals < 0] == 0
            
            if depths is None:
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0

            # skip sites where any of the queried individuals have less than `min_qual`
            if any([quals[i] < min_qual for i in (ki, si, mi, di, smi, sdi)]): continue

            # skip sites where any of the queried individuals have less than `min_depth` 
            if any([(ref_depths[i] + alt_depths[i]) < min_depth for i in (ki, si, mi, di, smi, sdi)]): 
                continue

            # grandparents must have evidence of putative variants.
            # this is mostly to speed up analysis, next chunk of code looks for
            # high quality evidence
            if sum([1 if gt in (HET, HOM_ALT) else 0 for gt in (m_gt, d_gt, sm_gt, sd_gt)]) < 1:
                continue
            
            grandparents_with_v = []

            # check grandparents for **HQ** evidence of the variant
            for gp in (kid.mom, kid.dad, spouse.mom, spouse.dad):
                idx = sample_dict[gp.sample_id]
                gt = (gts[idx][0], gts[idx][1])
                gp_qual = quals[idx]
                gp_rd, gp_ad = ref_depths[idx], alt_depths[idx]
                if gp_rd + gp_ad < min_depth: continue
                if gp_qual < min_qual: continue
                if gt in (HET, HOM_ALT): 
                    grandparents_with_v.append(gp.sample_id)
                else: continue
            
            if len(grandparents_with_v) == 0: continue

            family_idxs = [sample_dict[k.sample_id] for k in all_samples if k.family_id == kid.family_id]
            gt_to_filter = HET
            # NOTE: what about if the X variant is HOM_ALT in the guy, it could still be present
            # as a DNM in a female in the cohort as a HET?
            if v.CHROM == 'X' and kid.sex == 'male': gt_to_filter = HOM_ALT

            possible_carriers = 0 # genotype only
            likely_carriers = 0 # must meet QUAL and AB threshold

            for sample in all_samples:
                # ignore potential carriers in the same family as the proband
                if sample.family_id == kid.family_id: continue
                sample_idx = sample_dict[sample.sample_id]
                genotype = gts[sample_idx]
                genotype_ = (genotype[0], genotype[1])
                if genotype_ != gt_to_filter: continue
                possible_carriers += 1
                sample_ref = ref_depths[sample_idx]
                sample_alt = alt_depths[sample_idx]
                sample_total = float(sample_ref + sample_alt)
                sample_qual = quals[sample_idx]

                if sample_qual < 10: continue
                if sample_total < 10: continue
                if gt_to_filter == HOM_ALT:
                    if sample_alt / sample_total < 0.75: continue
                elif gt_to_filter == HET:
                    if not 0.25 <= (sample_alt / sample_total) <= 0.75: continue
                likely_carriers += 1
            possible_carriers -= likely_carriers
                
            # check F2s for HQ evidence of the variant
            f2s_with_v = []
            f2s_with_v_dp = []
            f2s_with_v_ab = []
            f2s_with_v_gt = []
            f2s_with_v_gq = []
            for child in children:
                ci = sample_dict[child.sample_id]
                gt = (gts[ci][0], gts[ci][1])
                rd, ad, qual, ab = ref_depths[ci], alt_depths[ci], quals[ci], get_ab(ref_depths, alt_depths, ci)
                if rd + ad < min_depth: continue
                if qual < min_qual: continue
                if gt not in (HET, HOM_ALT): continue
                f2s_with_v.append(child)
                f2s_with_v_dp.append(rd + ad)
                f2s_with_v_ab.append(ab)
                f2s_with_v_gt.append(gt)
                f2s_with_v_gq.append(qual)
            
            if len(f2s_with_v) == 0: continue

            ret.extend(variant_info(v, kid, sample_dict, possible_carriers=possible_carriers,
                                                         likely_carriers=likely_carriers,
                                                         grandparents_with_v=grandparents_with_v, 
                                                         f2s_with_v=f2s_with_v,
                                                         f2s_with_v_dp=f2s_with_v_dp,
                                                         f2s_with_v_ab=f2s_with_v_ab,
                                                         f2s_with_v_gt=f2s_with_v_gt,
                                                         f2s_with_v_gq=f2s_with_v_gq,
                                                         f1_gt=k_gt, f1_spouse_gt=s_gt,
                                                         f1_mom_gt=m_gt, f1_dad_gt=d_gt,
                                                         f1_spouse_mom_gt=sm_gt,
                                                         f1_spouse_dad_gt=sd_gt,
                                                         alt_i=alt_index))
            
    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    #if len(ret) == 1: return ret[0]
    if len(ret) > 0: 
        return ret

def variant_info(v, kid, sample_dict,
                     f1_gt=None, f1_spouse_gt=None,
                     f1_mom_gt=None, f1_dad_gt=None,
                     f1_spouse_mom_gt=None,
                     f1_spouse_dad_gt=None,
                     grandparents_with_v=[],
                     f2s_with_v_dp=[],
                     f2s_with_v_ab=[],
                     f2s_with_v_gt=[],
                     f2s_with_v_gq=[],
                     f2s_with_v=[],
                     interval_type=None, alt_i=None,
                     possible_carriers=None, likely_carriers=None):

    gt_dict = {(0, 0): 0, (0, 1): 1, (0, 2): 1, (1, 1): 2, (2, 2): 2}
    
    ac = v.INFO.get('AC')
    if type(ac) is tuple:
        ac = ac[alt_i]

    
    yield OrderedDict((
        ("chrom", v.CHROM),
        ("start", v.start),
        ("end", v.end),
        ("sample_id", kid.sample_id),
        ("family_id", kid.family_id),
        ("paternal_id", kid.paternal_id),
        ("maternal_id", kid.maternal_id),
        ("sample_sex", kid.sex),
        ("ref", v.REF),
        ("alt", v.ALT[alt_i]),
        ("alt_i", (str(alt_i + 1) + '/' + str(len(v.ALT)))), 
        ("ac", ac),
        ("an", v.INFO.get('AN')),
        ("possible_carriers", possible_carriers),
        ("likely_carriers", likely_carriers),
        ("mq_vcf", v.INFO.get("MQ")),
        ("mapq_ranksum_vcf", v.INFO.get('MQRankSum')),
        ("qual_by_depth", v.INFO.get('QD')),
        ("num_hom_alt", v.num_hom_alt),
        ("filter", v.FILTER or "PASS"),
        ("f1_gt", gt_dict[f1_gt]),
        ("f1_spouse_gt", gt_dict[f1_spouse_gt]),
        ("f1_mom_gt", gt_dict[f1_mom_gt]),
        ("f1_dad_gt", gt_dict[f1_dad_gt]),
        ("f1_spouse_mom_gt", gt_dict[f1_spouse_mom_gt]),
        ("f1_spouse_dad_gt", gt_dict[f1_spouse_dad_gt]),
        ("f2s_with_v", ','.join([x.sample_id for x in f2s_with_v])),
        ("f2s_with_v_dp", ','.join([str(x) for x in f2s_with_v_dp])),
        ("f2s_with_v_ab", ','.join([str(x) for x in f2s_with_v_ab])),
        ("f2s_with_v_gt", ','.join([str(gt_dict[x]) for x in f2s_with_v_gt])),
        ("f2s_with_v_gq", ','.join([str(x) for x in f2s_with_v_gq])),
        ("grandparents_with_v", ','.join(grandparents_with_v)),
        ("qual_mean", "%.1f" % np.mean(v.gt_quals)),
        ("call_rate", "%.3f" % v.call_rate)
        ))

def denovo(v, sample_dict, p0, f1s, f2s, max_alts_in_parents=2,
        min_depth=1, 
        min_depth_percentile=3,
        min_qual=1,
        multigen=False,
        _use_cohort_filters=True,
        exclude=None,
        simons_prefix=False,
        HET=1):
    if not variant_prefilter(v, 1): return None
    return get_denovo(v, sample_dict, p0, f1s, f2s,
            min_depth=min_depth,
            min_qual=min_qual,
            multigen=multigen,
            exclude=exclude,
            _use_cohort_filters=_use_cohort_filters,
            simons_prefix=simons_prefix,
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

    simons_prefix = False
    if args.simons_prefix:
        simons_prefix = True
    n_dn = 0
    for i, v in enumerate(vcf(str(args.chrom)) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)

        d = denovo(v, sample_dict, p0, f1s, f2s,
                exclude=exclude, _use_cohort_filters=_use_cohort_filters,
                simons_prefix=simons_prefix, min_qual=args.min_qual, min_depth=args.min_depth)
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
    p.add_argument("-simons_prefix", help="require all sample names to begin with the 'SSC' prefix", action="store_true")
    p.add_argument("-min_qual", help="minimum genotype quality to filter on [10]", default=10, type=int)
    p.add_argument("-min_depth", default=10, type=int)
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
