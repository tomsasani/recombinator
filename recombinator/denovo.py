from __future__ import print_function, absolute_import, division
import os
import sys
import time
from collections import OrderedDict

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

    if len(v.REF) > 4: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if len(v.ALT[0]) > 6: return False
    if v.FILTER is not None and not tranche99(v.FILTER) : return False

    if v.QUAL < min_variant_qual: return False
    return True

def get_denovo(v, samples, f1s, f2s, max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_qual=1,
        exclude=None,
        _use_cohort_filters=True,
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

    if _use_cohort_filters and v.num_het > 2: return None
    if _use_cohort_filters and v.num_hom_alt > 1: return None

    ret = []
    gts = v.gt_types

    depths = None
    ref_depths = None

    # loop over alternate alleles
    for k in range(1, len(v.ALT) + 1):
        alt_depths = None

        # and then loop over the f1s. 
        for f1 in f1s:
            f1_idx = samples[f1.sample_id]
            if gts[f1_idx] != HET: continue
            # check that parents are hom-ref
            if f1.mom.sample_id == '8477' or f1.dad.sample_id == '8477':
                continue
            mi, di = samples[f1.mom.sample_id], samples[f1.dad.sample_id]
            # Added try/except since sample '8477' apparently doesn't exist in a pedigree
            for pi in (mi, di):
                if not (gts[pi] == 0 or gts[pi] == 2): continue
            # List of peddy.Samples() that represent children of the f1 we're looping over.
            if f2s is not None:
                f2s = [k for k in f2s if f1.sample_id in [k.mom.sample_id, k.dad.sample_id]]
                f2_idxs = [samples[k.sample_id] for k in f2s]

            if depths is None:
                #depths = v.format('AD', int)
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0

            if alt_depths[[mi, di]].sum() > max_alts_in_parents: continue

            f1_alt = alt_depths[f1_idx]
            f1_ref = ref_depths[f1_idx]

            # depth filtering.
            if f1_alt + f1_ref < min_depth + 1: continue
            if f1_alt < min_depth / 2: continue
            if ref_depths[di] + alt_depths[di] < min_depth - 1: continue
            if ref_depths[mi] + alt_depths[mi] < min_depth - 1: continue

            if np.mean(alt_depths + ref_depths) > max_mean_depth:
                continue

            # if there are too many alts outside this kid. skip
            # however, ignore the potential alts (i.e., transmitted DNM) in the f1's
            # children
            alt_sum = sum([alt_depths[i] for i in range(len(alt_depths)) if i not in f2_idxs]) - f1_alt
            ref_sum = sum([ref_depths[i] for i in range(len(ref_depths)) if i not in f2_idxs]) - f1_ref

            # this check is less stringent than the p-value checks below but
            # avoids some compute.
            if _use_cohort_filters and alt_sum > len(samples) * 0.01:
                continue

            # balance evidence in kid and alts in entire cohort.
            if _use_cohort_filters and alt_sum >= f1_alt: continue

            if _use_cohort_filters and alt_sum > 0:
                if f1_alt < 6: continue

            # via Tom Sasani.
            palt = ss.binom_test([alt_sum, ref_sum], 
                    p=0.0002,
                    alternative="greater")
            if _use_cohort_filters and palt < min_allele_balance_p: continue

            pab = ss.binom_test([f1_ref, f1_alt])
            if pab < min_allele_balance_p: continue

            # TODO: check why some quals are 0 and if filtering on this improve
            #quals = v.gt_quals
            # accuracy.
            #quals[quals < 0] == 0
            #if quals[ki] < 1 or quals[mi] < 1 or quals[di] < 1: continue

            # stricter settings with FILTER
            if v.FILTER is not None:
                # fewer than 1 alt per 500-hundred samples.
                if _use_cohort_filters and alt_sum > 0.002 * len(samples): continue
                # no alts in either parent.
                if alt_depths[[mi, di]].sum() > 0: continue

            # at this point, the variant is a solid putative DNM. now need
            # to validate that it's in at least one of the f2s. all of the
            # previous filters were on a population basis, so for now we'll only
            # check if one of the kids is a HET
            f2_count = 0
            for f2 in f2s:
                mi, di = samples[f2.mom.sample_id], samples[f2.dad.sample_id]
                # Make sure that the F1's partner doesn't also have the same DNM
                if gts[mi] + gts[di] != 1: continue
                f2_idx = samples[f2.sample_id]
                if gts[f2_idx] == HET:
                    f2_count += 1
            if f2_count > 0:
                ret.extend(variant_info(v, f1, samples, pab, palt, alt_i=k))

    if exclude is not None and 0 != len(exclude[v.CHROM].search(v.start, v.end)):
        return None

    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    if len(ret) == 1: return ret[0]

def variant_info(v, kid, samples, pab=None, palt=None, alt_i=None):

    quals = v.gt_quals
    ki, mi, di = samples[kid.sample_id], samples[kid.mom.sample_id], samples[kid.dad.sample_id]
    depths = v.format('AD', int)
    ref_depths = depths[:, 0]
    all_alts = depths[:, 1:]
    for k in range(all_alts.shape[1]):
        if alt_i is not None and alt_i != (k + 1): continue
        alt_depths = all_alts[:, k]

        #ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals
        kid_ref, kid_alt = ref_depths[ki], alt_depths[ki]
        alt_sum = alt_depths.sum() - kid_alt
        if pab is None:
            pab = ss.binom_test([kid_ref, kid_alt])
        if palt is None:
            palt = ss.binom_test([alt_sum, ref_depths.sum() - kid_ref], p=0.0002,
                    alternative="greater")

        yield OrderedDict((
            ("chrom", v.CHROM),
            ("start", v.start),
            ("end", v.end),
            ("sample_id", f1.sample_id),
            ("family_id", f1.family_id),
            ("ref", v.REF),
            ("alt", v.ALT[k]),
            ("alt_i", "%d/%d" % (k, len(v.ALT))),
            ("filter", v.FILTER or "PASS"),
            ("pab", "%.3g" % pab),
            ("palt", "%.3g" % palt),
            ("paternal_id", f1.paternal_id),
            ("maternal_id", f1.maternal_id),
            ("kid_ref_depth", f1_ref),
            ("kid_alt_depth", f1_alt),
            ("kid_total_depth", f1_ref + f1_alt),
            ("mom_ref_depth", ref_depths[mi]),
            ("mom_alt_depth", alt_depths[mi]),
            ("mom_total_depth", ref_depths[mi] + alt_depths[mi]),
            ("dad_ref_depth", ref_depths[di]),
            ("dad_alt_depth", alt_depths[di]),
            ("dad_total_depth", ref_depths[di] + alt_depths[di]),
            ("kid_qual", quals[f1_idx]),
            ("mom_qual", quals[mi]),
            ("dad_qual", quals[di]),
            #("p%d_depth" % min_depth_percentile, p5),
            ("depth_mean", "%.1f" % np.mean(np.array(ref_sum) + np.array(alt_sum))),
            ("depth_std", "%.1f" % np.std(np.array(ref_sum) + np.array(alt_sum))),
            ("qual_mean", "%.1f" % np.mean(quals)),
            ("call_rate", "%.3f" % v.call_rate),
            ("cohort_alt_depth", alt_sum),
            ))


def denovo(v, samples, f1s, f2s,  max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_depth_percentile=3,
        min_qual=1,
        _use_cohort_filters=True,
        exclude=None,
        HET=1):
    if not variant_prefilter(v, 10): return None
    return get_denovo(v, samples, f1s, f2s, max_alts_in_parents=max_alts_in_parents,
            min_depth=min_depth,
            max_mean_depth=max_mean_depth,
            min_allele_balance_p=min_allele_balance_p,
            min_qual=min_qual,
            exclude=exclude,
            _use_cohort_filters=_use_cohort_filters,
            HET=HET)


def write_denovo(d, fh, _line={}):
    if _line.get(fh.name, 0) == 0:
        _line[fh.name] = 1
        print("#" + "\t".join(d.keys()), file=fh)
    print("\t".join(map(str, d.values())), file=fh)

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

    samples_lookup = {v: i for i, v in enumerate(vcf.samples)}
    kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    p0 = [k for k in ped.samples() if k.mom is None and k.dad is None]
    f1s = [k for k in kids if k.mom.mom is None and k.dad.dad is None]
    f2s = [k for k in kids if k not in f1s]

    _use_cohort_filters = False if os.environ.get("DN_NO_COHORT") else True

    n_dn = 0
    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)

        d = denovo(v, samples_lookup, f1s, f2s, max_alts_in_parents=1,
                exclude=exclude, _use_cohort_filters=_use_cohort_filters)
        if d is not None:
            write_denovo(d, sys.stdout)
            n_dn += 1

    sys.stdout.flush()
    print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)


def main(argv):
    from argparse import ArgumentParser
    p = ArgumentParser()
    p.add_argument("--chrom", help="extract only this chromosome", default=None)
    p.add_argument("--exclude", help="regions to exclude (e.g. LCRs)")
    p.add_argument("--min-mean-depth", default=10, type=int)
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
