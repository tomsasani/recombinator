from __future__ import print_function, absolute_import, division
import os
import sys
import time
from collections import OrderedDict

import pysam
from peddy import Ped
import cyvcf2
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
    #if len(v.REF) > 10: return False
    #if len(v.ALT) > 2 or "*" in v.ALT: return False
    if "*" in v.ALT: return False
    #if len(v.ALT[0]) > 10: return False
    #if v.FILTER is not None and not tranche99(v.FILTER): return False

    if v.QUAL < min_variant_qual: return False
    return True

def validate_in_f2(v, sample_dict, gts, kid, f2s, f1_dnm_allele=1):
    """
    at this point, the variant is a solid putative DNM. now need
    to validate that it's in at least one of the f2s. all of the
    previous filters were on a population basis, so for now we'll only
    check if one of the kids is a HET
    """

    f1_is_mom = kid.sex == 'female'
    f1_is_dad = kid.sex == 'male'
    # define genotypes
    UNKNOWN = (-1, -1)
    HOM_REF = (0, 0)
    f1_dnm_allele=f1_dnm_allele[0] 
    quals = v.gt_quals
    f2_with_v = []
    for f2 in f2s:
        if kid.sample_id not in [f2.mom.sample_id, f2.dad.sample_id]: continue

        f2_idx = sample_dict[f2.sample_id]
        mi, di = sample_dict[f2.mom.sample_id], sample_dict[f2.dad.sample_id]
        
        m_gt, d_gt = (gts[mi][0], gts[mi][1]), (gts[di][0], gts[di][1])
        f2_gt = (gts[f2_idx][0], gts[f2_idx][1])
        if any([x == UNKNOWN for x in (m_gt, d_gt)]): continue

        # Make sure that the F1's partner doesn't also have the same DNM
        # Actually, this is possible? But extreeeeemely unlikely. This check only
        # works if we don't look at the X.
        if f1_is_mom and f1_dnm_allele in d_gt: continue
        elif f1_is_dad and f1_dnm_allele in m_gt: continue
        else: pass

        if v.CHROM != 'X':
            if list(f2_gt).count(f1_dnm_allele) != 1: continue
        # Special logic for inherited DNM on male X.
        elif v.CHROM == 'X':
            # Female F1 with DNM passes it on to male F2, must be
            # homozygous in the F2.
            if kid.sex == 'female' and f2.sex == 'male':
                if not all([x == f1_dnm_allele for x in list(f2_gt)]): continue 
            # Male F1 with DNM on X simply cannot pass it onto male F2.
            elif kid.sex == 'male' and f2.sex == 'male': continue
            else:
                if list(f2_gt).count(f1_dnm_allele) != 1: continue

        f2_with_v.append(f2)

    return f2_with_v

def get_denovo(v, sample_dict, p0, f1s, f2s,
        min_depth=5,
        max_mean_depth=400,
        min_qual=1,
        multigen=False,
        exclude=None,
        _use_cohort_filters=True,
        simons_prefix=False,
        eiee=False,
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

    ret = []
    gts = v.genotypes # NOTE: changes v.gt_types to v.genotypes for multiallelics
   
    all_samples = p0 + f1s + f2s
    all_children = f1s + f2s
    if multigen: kids = f1s
    else: kids = f2s
    if eiee or simons_prefix: kids = f1s
    if simons_prefix:
        all_children = [x for x in all_children if 'SSC' in x.sample_id]

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

            if v.CHROM == 'Y': continue # ignore Y for now
            # define genotypes
            UNKNOWN = (-1, -1)
            HOM_REF = (0, 0)
            
            k_gt, m_gt, d_gt = (gts[ki][0], gts[ki][1]),\
                               (gts[mi][0], gts[mi][1]),\
                               (gts[di][0], gts[di][1])

            # get a list of the parent genotypes, which we'll compare to the kid's genotype
            parent_gts = list(m_gt)
            parent_gts.extend(list(d_gt))
            
            # if the kid has a genotype that's not observed in either parent, it's de novo
            uniq_gt_in_kid = any([x not in parent_gts for x in list(k_gt)])
            f1_dnm_allele = [a for a in list(k_gt) if a not in parent_gts]

            if not uniq_gt_in_kid: continue
            if alt_index + 1 != f1_dnm_allele[0]: continue

            # can't be confident in any genotypes that are UNKNOWN
            if any([gt == UNKNOWN for gt in (k_gt, m_gt, d_gt)]): continue

            kid_is_hom = len(set(k_gt)) == 1
            kid_is_het = len(set(k_gt)) == 2

            # check that kid is a HET (or HOM_ALT on hemizygous chroms in males)
            if kid.sex == 'male' and v.CHROM == 'X':
                if not kid_is_hom: continue
            else:
                if not kid_is_het: continue

            quals = v.gt_quals
            quals[quals < 0] == 0
            if quals[ki] < min_qual or quals[mi] < min_qual or quals[di] < min_qual: continue
           
            if depths is None:
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0

            m_alt, m_ref = alt_depths[mi], ref_depths[mi]
            d_alt, d_ref = alt_depths[di], ref_depths[di]
            k_alt, k_ref = alt_depths[ki], ref_depths[ki]
            
            # check total depth in mom, dad, and kid, and for lack of ALT support
            if v.CHROM == 'X': norm_depth_by = 2.
            else: norm_depth_by = 1.
            if (m_ref + m_alt) < (min_depth / norm_depth_by): continue
            if (d_ref + d_alt) < (min_depth / norm_depth_by): continue 
            if (k_alt + k_ref) < (min_depth / norm_depth_by): continue

            ab = k_alt / float(k_alt + k_ref)

            # if we're looking in a third generation, kids (F1s) have children
            children = []
            if multigen:
                children = [k for k in f2s if kid.sample_id in (k.dad.sample_id, k.mom.sample_id)]
                       
            family_idxs = [sample_dict[k.sample_id] for k in all_samples if k.family_id == kid.family_id]
            #print ('allele to filter is {}'.format(allele_to_filter))
            # NOTE: what about if the X variant is HOM_ALT in the guy, it could still be present
            # as a DNM in a female in the cohort as a HET?
            #if v.CHROM == 'X' and kid.sex == 'male': gt_to_filter = k_gt

            possible_carriers = 0 # genotype only
            likely_carriers = 0 # must meet QUAL and AB threshold

            allele_to_filter = f1_dnm_allele[0]

            for sample in all_samples:
                if sample.family_id == kid.family_id: continue
                sample_idx = sample_dict[sample.sample_id]
                genotype = gts[sample_idx]
                genotype_ = (genotype[0], genotype[1])
                if allele_to_filter not in genotype_: continue
                possible_carriers += 1
                sample_ref = ref_depths[sample_idx]
                sample_alt = alt_depths[sample_idx]
                sample_total = float(sample_ref + sample_alt)
                sample_qual = quals[sample_idx]

                if sample_qual < 20: continue
                if sample_total < 10: continue
                if kid_is_hom:
                    if sample_alt / sample_total < 0.75: continue
                elif kid_is_het:
                    if not 0.25 <= (sample_alt / sample_total) <= 0.75: continue
                likely_carriers += 1
            possible_carriers -= likely_carriers
            
            # check for inheritance of the DNM in a third generation
            putative_f2_with_v = []
            if multigen:
                putative_f2_with_v = validate_in_f2(v, sample_dict, gts, kid, children, f1_dnm_allele=f1_dnm_allele)

            f2_with_v = []
            f2_with_v_ab = []
            f2_with_v_gq = []
            f2_with_v_dp = []
            # impose lenient filters on individuals that inherit the DNM
            # track allele balance in transmission events to look for mosaics
            for f2 in putative_f2_with_v:
                f2_idx = sample_dict[f2.sample_id]
                total_depth = alt_depths[f2_idx] + ref_depths[f2_idx]
                ab = alt_depths[f2_idx] / float(total_depth)
                if ab == 0.: continue
                f2_with_v.append(f2.sample_id)
                f2_with_v_ab.append(ab)
                f2_with_v_gq.append(quals[f2_idx])
                f2_with_v_dp.append(total_depth)

            ret.extend(variant_info(v, kid, sample_dict, f2_with_v=f2_with_v, 
                                                         f2_with_v_ab=f2_with_v_ab,
                                                         f2_with_v_gq=f2_with_v_gq,
                                                         f2_with_v_dp=f2_with_v_dp,
                                                         alt_i=alt_index,
                                                         possible_carriers=possible_carriers,
                                                         likely_carriers=likely_carriers))
    # assuming you have multiallelics        
    if len(ret) > 0: 
        return ret

def variant_info(v, kid, sample_dict, f2_with_v=[], f2_with_v_ab=[], f2_with_v_gq=[], f2_with_v_dp=[],
                    interval_type=None, alt_i=None, possible_carriers=0, likely_carriers=0):

    ki, mi, di = sample_dict[kid.sample_id], sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

    f2_with_v_idx = [sample_dict[f2] for f2 in f2_with_v]

    ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals

    k_ref, k_alt = ref_depths[ki], alt_depths[ki]

    alt_sum = sum([alt_depths[i] for i in range(len(alt_depths)) if i not in f2_with_v_idx]) - k_alt
    ref_sum = sum([ref_depths[i] for i in range(len(ref_depths)) if i not in f2_with_v_idx]) - k_ref

    ac = v.INFO.get('AC')
    mleaf_vcf = v.INFO.get('MLEAF')
    mleac_vcf = v.INFO.get('MLEAC')
    if type(ac) is tuple:
        ac = ac[alt_i]
        mleaf_vcf = mleaf_vcf[alt_i]
        mleac_vcf = mleac_vcf[alt_i]


    f2s, f2_abs, f2_dps, f2_gqs = 'None', 'None', 'None', 'None'

    if len(f2_with_v) > 0:
        f2s = ','.join([str(x) for x in f2_with_v])
        f2_abs = ','.join([str(x) for x in f2_with_v_ab])
        f2_dps = ','.join([str(x) for x in f2_with_v_dp])
        f2_gqs = ','.join([str(x) for x in f2_with_v_gq])


    yield OrderedDict((
        ("chrom", v.CHROM),
        ("start", v.start),
        ("end", v.end),
        ("sample_id", kid.sample_id),
        ("family_id", kid.family_id),
        ("sample_sex", kid.sex),
        ("ref", v.REF),
        ("alt", v.ALT[alt_i]),
        ("alt_i",  (str(alt_i + 1) + '/' + str(len(v.ALT)))),
        ("ac", ac),
        ("an", v.INFO.get('AN')),
        ("possible_carriers", possible_carriers),
        ("likely_carriers", likely_carriers),
        ("num_hom_alt", v.num_hom_alt),
        ("mq0_vcf", v.INFO.get("MQ0")),
        ("mq_vcf", v.INFO.get("MQ")),
        ("mapq_ranksum_vcf", v.INFO.get('MQRankSum')),
        ("base_q_ranksum_vcf", v.INFO.get("BaseQRankSum")),
        ("clipping_ranksum_vcf", v.INFO.get("ClippingRankSum")),
        ("excess_het_vcf", v.INFO.get("ExcessHet")),
        ("strand_bias_vcf", v.INFO.get("FS")),
        ("inbreeding_coef_vcf", v.INFO.get("InbreedingCoeff")),
        ("mleac_vcf", mleac_vcf),
        ("mleaf_vcf", mleaf_vcf),
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
        ("cohort_alt_depth", alt_sum),
        ("f2_with_variant", f2s),
        ("f2_with_variant_ab", f2_abs),
        ("f2_with_variant_gq", f2_gqs),
        ("f2_with_variant_dp", f2_dps)
        ))

def denovo(v, sample_dict, p0, f1s, f2s,
        min_depth=10, 
        max_mean_depth=400,
        min_depth_percentile=3,
        min_qual=1,
        multigen=False,
        _use_cohort_filters=True,
        exclude=None,
        simons_prefix=False,
        eiee=False,
        HET=1):
    if not variant_prefilter(v, 1): return None
    return get_denovo(v, sample_dict, p0, f1s, f2s,
            min_depth=min_depth,
            max_mean_depth=max_mean_depth,
            min_qual=min_qual,
            multigen=multigen,
            exclude=exclude,
            _use_cohort_filters=_use_cohort_filters,
            simons_prefix=simons_prefix,
            eiee=eiee,
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
        if i % 10000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)

        d = denovo(v, sample_dict, p0, f1s, f2s, min_depth=args.min_depth,
                min_qual=args.min_qual,
                multigen=args.multigen, eiee=args.eiee,
                exclude=exclude, _use_cohort_filters=_use_cohort_filters,
                simons_prefix=simons_prefix)
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
    p.add_argument("-min_depth", help="min depth in kid and both parents [5]", type=int, default=5)
    p.add_argument("-min_qual", help="minimum GQ in kid [1]", type=int, default=1)
    p.add_argument("-simons_prefix", help="require all sample names to begin with the 'SSC' prefix", action="store_true")
    p.add_argument("-eiee", help="placeholder for eiee calling", action="store_true")
    p.add_argument("-multigen", help="use a third generation of samples to validate DNMs", action="store_true")
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
