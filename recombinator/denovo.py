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
    # Current hard filter on indel mutations.
    if len(v.REF) > 1: return False
    if len(v.ALT) > 1 or "*" in v.ALT: return False
    # Current hard filter on indel mutations.
    if len(v.ALT[0]) > 1: return False
    if v.FILTER is not None and not tranche99(v.FILTER) : return False

    if v.QUAL < min_variant_qual: return False
    return True

def validate_in_f2(v, sample_dict, gts, kid, f2s):
    """
    at this point, the variant is a solid putative DNM. now need
    to validate that it's in at least one of the f2s. all of the
    previous filters were on a population basis, so for now we'll only
    check if one of the kids is a HET
    """
    f2_with_v = []
    # Calculate the naive probability that we wouldn't observe this 
    # candidate DNM in one of the F2s.
    p_not_seen = 0.5 ** (len(f2s))
    for f2 in f2s:
        mi, di = sample_dict[f2.mom.sample_id], sample_dict[f2.dad.sample_id]
        assert kid.sample_id in [f2.mom.sample_id, f2.dad.sample_id]
        # Make sure that the F1's partner doesn't also have the same DNM
        # Actually, this is possible? But extreeeeemely unlikely. This check only
        # works if we don't look at the X.
        if v.CHROM != 'X' and (gts[mi] + gts[di]) > 1:
            continue
        f2_idx = sample_dict[f2.sample_id]
        # Special logic for inherited DNM on male X.
        if v.CHROM == 'X':
            # Female F1 with DNM passes it on to male F2, must be
            # homozygous in the F2.
            if kid.sex == 'female' and f2.sex == 'male':
                if gts[f2_idx] == 3:
                    f2_with_v.append(f2)
            # Male F1 with DNM on X simply cannot pass it onto male F2.
            elif kid.sex == 'male' and f2.sex == 'male':
                continue
            else:
                if gts[f2_idx] == 1:
                    f2_with_v.append(f2)
        else:
            if gts[f2_idx] == 1:
                f2_with_v.append(f2)

    return f2_with_v

def enriched_for_splitters(v, f1, split_limit=5,
    prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):

    split_count = 0

    bam = pysam.AlignmentFile(prefix + f1.sample_id + ".bam", "rb")
    reads = []
    for read in bam.fetch(v.CHROM, int(v.POS) - 15, int(v.POS) + 15):
        try:
            if read.query_name not in reads and len(read.get_tag("SA")) > 0:
                reads.append(read.query_name)
                split_count += 1
                continue
        except KeyError:
            continue

    if split_count >= split_limit:
        return True

def eval_alignments(v, f1, 
    prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):

    m_bam = pysam.AlignmentFile(prefix + f1.mom.sample_id + ".bam", "rb")
    m_ab = alignment_ab(v, m_bam)
    d_bam = pysam.AlignmentFile(prefix + f1.dad.sample_id + ".bam", "rb")
    d_ab = alignment_ab(v, d_bam)

    return (m_ab, d_ab)

def alignment_ab(v, bam):
    ref_sum, alt_sum = 0, 0
    ref = v.REF
    if isinstance(v.ALT, list):
        alt = v.ALT[0]
    else: alt = v.ALT
    for pileupcolumn in bam.pileup(v.CHROM, int(v.POS) - 1, int(v.POS) + 1):
        if pileupcolumn.pos != int(v.POS) - 1:
            continue
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                q_base = pileupread.alignment.query_sequence[pileupread.query_position]
                if q_base.upper() == ref.upper():
                    ref_sum += 1
                elif q_base.upper() == alt.upper():
                    alt_sum += 1
                else: continue
    ab = alt_sum / float(ref_sum + alt_sum) 
    return ab

def enriched_for_bad_mapping(v, f1, mapq0_lim=5,
    prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):
    bam = pysam.AlignmentFile(prefix + f1.sample_id + ".bam", "rb")

    mapq0 = 0
    for read in bam.fetch(v.CHROM, int(v.POS) - 10, int(v.POS) + 10):
        if read.mapping_quality == 0:
            mapq0 += 1
    if mapq0 >= mapq0_lim:
        return True

def get_denovo(v, sample_dict, kids, f2s=None,
        max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_cohort_enrich_p=0.01,
        min_qual=1,
        min_p_ab=0.1,
        exclude=None,
        one_sample=None,
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

    if exclude is not None:
        if len(exclude[v.CHROM].search(v.start, v.end)) > 0:
            return None

    #if _use_cohort_filters and v.num_hom_alt > 1: return None

    ret = []
    gts = v.gt_types

    depths = None
    ref_depths = None

    # loop over alternate alleles
    for k in range(1, len(v.ALT) + 1):
        alt_depths = None

        # and then loop over the f1s. 
        for kid in kids:
            k_idx = sample_dict[kid.sample_id]

            # If we're looking for DNMs in the F2 generation,
            # we don't expect to see this DNM *anywhere* else in the cohort.
            if f2s is None and v.num_het > 1:
                continue

            # Special logic for X chromosome.
            if kid.sex == 'male' and v.CHROM == 'X':
                if gts[k_idx] != 3: continue
            else:
                if gts[k_idx] != 1: continue
            # Added check since sample '8477' apparently doesn't exist in our ceph pedigree
            if kid.mom.sample_id == '8477' or kid.dad.sample_id == '8477':
                continue
            mi, di = sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]
            # check that parents are hom-ref
            if not (gts[di] == 0 and gts[mi] == 0):
                continue

            ac_in_f2s = 0
            f2_idxs = []

            # List of peddy.Samples() that represent children of the f1 we're looping over.
            # Assuming there are any, of course.
            if f2s is not None:
                f2s = [k for k in f2s if kid.sample_id in [k.mom.sample_id, k.dad.sample_id]]
                f2_idxs = [sample_dict[k.sample_id] for k in f2s]

                # Total number of ALT alleles in F2s.
                ac_dict = {3:2, 2:0, 1:1, 0:0}
                ac_in_f2s = sum([ac_dict[gts[i]] for i in f2_idxs])

            if depths is None:
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0

            if alt_depths[[mi, di]].sum() > max_alts_in_parents: 
                continue

            if ref_depths[mi] + alt_depths[mi] < min_depth:
                continue
            if ref_depths[di] + alt_depths[di] < min_depth:
                continue

            k_alt, k_ref = alt_depths[k_idx], ref_depths[k_idx]
            if k_alt + k_ref < min_depth:
                continue
            if np.mean(alt_depths + ref_depths) > max_mean_depth:
                continue

            # if there are too many alts outside this kid, skip. however,
            # ignore the potential alts (i.e., transmitted DNM) in the f1's children
            alt_sum = sum([alt_depths[i] for i in range(len(alt_depths)) if i not in f2_idxs]) - k_alt
            ref_sum = sum([ref_depths[i] for i in range(len(ref_depths)) if i not in f2_idxs]) - k_ref
            # this check is less stringent than the p-value checks below but
            # avoids some compute.
            # NOTE: got rid of this
            #if _use_cohort_filters and alt_sum > len(samples) * 0.01: continue

            # balance evidence in kid and alts in entire cohort.
            if _use_cohort_filters and alt_sum >= k_alt: continue

            # NOTE: got rid of below filter
            #if _use_cohort_filters and alt_sum > 0:
            #    if f1_alt < 6: continue

            # via Tom Sasani.
            # NOTE: changed this to do the test on number of alt ALLELEs, not depth
            # calculate enrichment of ALT alleles in the cohort. [x, y] where x is 
            # the number of ALT alleles in the cohort (outside the F1 we're currently looking at and
            # the F1's kids), and y is the total number of alleles in the cohort, outside 
            # of the F2s and the F1 we're looking at

            #alt_in_cohort = int(v.INFO.get('AC')) - 1 - ac_in_f2s
            #total_alleles_in_cohort = int(v.INFO.get('AN')) - (len(f2s) * 2) - 2

            palt = ss.binom_test([alt_sum, ref_sum],  
                    p=0.002, # expected probability of observing an alt
                               # read in the population
                    alternative="greater") # one-sided
            if _use_cohort_filters and palt < min_cohort_enrich_p: continue

            pab = ss.binom_test([k_ref, k_alt], p=0.5, alternative='two-sided') 
            if pab < min_allele_balance_p: continue

            # TODO: check why some quals are 0 and if filtering on this improve
            # accuracy.
            quals = v.gt_quals
            quals[quals < 0] == 0
            if quals[k_idx] < 1 or quals[mi] < 1 or quals[di] < 1: continue

            # stricter settings with FILTER
            # Note: seems as though FILTER is always None...
            if v.FILTER is not None:
                # fewer than 1 alt per 500-hundred samples.
                #if _use_cohort_filters and alt_sum > 0.002 * len(samples): continue
                if _use_cohort_filters and int(v.INFO.get('AC')) > 0.002 * len(sample_dict): continue

            # Confirm that parents don't have excess of ALT reads in the BAM.
            m_ab, d_ab = eval_alignments(v, kid)
            if (m_ab > min_p_ab or d_ab > min_p_ab):
                continue
            # Check for enrichment of split alignments surrounding the putative DNM.
            # If True, indicative of a problematic region for alignments.
            if enriched_for_splitters(v, kid, split_limit=5):
                continue
            # Check for enrichment of 0-quality alignments, similar to check above.
            if enriched_for_bad_mapping(v, kid, mapq0_lim=5):
                continue

            if f2s is not None:
                f2_with_v = validate_in_f2(v, sample_dict, gts, kid, f2s)
            else: f2_with_v = []
            p_dnm_in_f2s = 1 - 0.5 ** (len(f2s) - len(f2_with_v))

            ret.extend(variant_info(v, kid, sample_dict, f2_with_v=f2_with_v, 
                                                        pab=pab, palt=palt, alt_i=k,
                                                        p_dnm_in_f2s=p_dnm_in_f2s))


    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    if len(ret) == 1: return ret[0]

def variant_info(v, kid, sample_dict, f2_with_v=[], 
                    interval_type=None, pab=None, palt=None, alt_i=None, p_dnm_in_f2s=None):

    k_idx, mi, di = sample_dict[kid.sample_id], sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

    f2_with_v_idx = [sample_dict[f2.sample_id] for f2 in f2_with_v]
    #depths = v.format('AD', int)
    #ref_depths = depths[:, 0]
    #all_alts = depths[:, 1:]

    ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals

    k_ref, k_alt = ref_depths[k_idx], alt_depths[k_idx]

    alt_sum = sum([alt_depths[i] for i in range(len(alt_depths)) if i not in f2_with_v_idx]) - k_alt
    ref_sum = sum([ref_depths[i] for i in range(len(ref_depths)) if i not in f2_with_v_idx]) - k_ref

    #for k in range(all_alts.shape[1]):
    
    for k in range(len(v.ALT)):

        # Below line prevents DNM calling, need to figure out why.
        #if alt_i is not None and alt_i != (k + 1): continue
        #alt_depths = all_alts[:, k]
    
        yield OrderedDict((
            ("chrom", v.CHROM),
            ("start", v.start),
            ("end", v.end),
            ("sample_id", kid.sample_id),
            ("family_id", kid.family_id),
            ("sample_sex", kid.sex),
            ("ref", v.REF),
            ("alt", v.ALT[k]),
            ("alt_i", "%d/%d" % (k, len(v.ALT))),
            ("filter", v.FILTER or "PASS"),
            ("pab", "%.3g" % pab),
            ("palt", "%.3g" % palt),
            ("paternal_id", kid.paternal_id),
            ("maternal_id", kid.maternal_id),
            ("kid_ref_depth", k_ref),
            ("kid_alt_depth", k_alt),
            ("kid_total_depth", k_ref + k_alt),
            ("mom_ref_depth", ref_depths[mi]),
            ("mom_alt_depth", alt_depths[mi]),
            ("mom_total_depth", ref_depths[mi] + alt_depths[mi]),
            ("dad_ref_depth", ref_depths[di]),
            ("dad_alt_depth", alt_depths[di]),
            ("dad_total_depth", ref_depths[di] + alt_depths[di]),
            ("kid_qual", quals[k_idx]),
            ("mom_qual", quals[mi]),
            ("dad_qual", quals[di]),
            #("p%d_depth" % min_depth_percentile, p5),
            ("depth_mean", "%.1f" % np.mean(ref_depths + alt_depths)),
            ("depth_std", "%.1f" % np.std(ref_depths + alt_depths)),
            ("qual_mean", "%.1f" % np.mean(quals)),
            ("call_rate", "%.3f" % v.call_rate),
            ("cohort_alt_depth", alt_sum),
            ("f2_with_variant", ','.join([str(x.sample_id) for x in f2_with_v])),
            ("p_dnm_in_f2s", p_dnm_in_f2s),
            ))


def denovo(v, sample_dict, kids, f2s=None, max_alts_in_parents=5,
        min_depth=10,
        max_mean_depth=400,
        min_allele_balance_p=0.05,
        min_cohort_enrich_p=0.01,
        min_depth_percentile=3,
        min_qual=1,
        _use_cohort_filters=True,
        exclude=None,
        HET=1):
    if not variant_prefilter(v, 10): return None
    return get_denovo(v, sample_dict, kids, f2s=f2s, max_alts_in_parents=max_alts_in_parents,
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

    sample_dict = {v: i for i, v in enumerate(vcf.samples)}
    kids = [k for k in ped.samples() if k.mom is not None and k.dad is not None]

    p0 = [k for k in ped.samples() if k.mom is None and k.dad is None]
    f1s = [k for k in kids if k.mom.mom is None and k.dad.dad is None and k.mom.dad is None and k.dad.mom is None]
    f2s = [k for k in kids if k not in f1s]

    # If we're looking for DNM in the F2s.
    if args.f2:
        kids = f2s
        f2s = None
    else: kids = f1s
    _use_cohort_filters = False if os.environ.get("DN_NO_COHORT") else True

    n_dn = 0
    for i, v in enumerate(vcf(args.chrom) if args.chrom else vcf, start=1):
        if i % 100000 == 0:
            print(" called %d de-novos out of %d variants" % (n_dn, i), file=sys.stderr)
        d = denovo(v, sample_dict, kids, f2s=f2s,
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
    p.add_argument("-f2", help="if we're looking for DNM in the F2s", action="store_true")
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
