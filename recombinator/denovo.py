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
    if len(v.REF) > 4: return False
    if len(v.ALT) > 2 or "*" in v.ALT: return False
    if len(v.ALT[0]) > 6: return False
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
    quals = v.gt_quals
    f2_with_v = []
    for f2 in f2s:
        f2_idx = sample_dict[f2.sample_id]
        mi, di = sample_dict[f2.mom.sample_id], sample_dict[f2.dad.sample_id]
        assert kid.sample_id in [f2.mom.sample_id, f2.dad.sample_id]
        # Make sure that the F1's partner doesn't also have the same DNM
        # Actually, this is possible? But extreeeeemely unlikely. This check only
        # works if we don't look at the X.
        if v.CHROM != 'X' and (gts[mi] + gts[di]) > 1: continue
        if v.CHROM == 'Y' and f2.sex == 'female': continue
        # Special logic for inherited DNM on male X.
        if v.CHROM == 'X':
            # Female F1 with DNM passes it on to male F2, must be
            # homozygous in the F2.
            if kid.sex == 'female' and f2.sex == 'male':
                if gts[f2_idx] == 3: f2_with_v.append(f2)
            # Male F1 with DNM on X simply cannot pass it onto male F2.
            elif kid.sex == 'male' and f2.sex == 'male': continue
            else:
                if gts[f2_idx] == 1: f2_with_v.append(f2)
        else:
            if gts[f2_idx] == 1: f2_with_v.append(f2)

    return f2_with_v

def query_alignments(v, sample_id, prefix="/scratch/ucgd/lustre/ugpuser/Repository/AnalysisData/2016/A414/16-08-06_WashU-Yandell-CEPH/UGP/Data/PolishedBams/"):

    from collections import namedtuple, defaultdict
    import scipy.stats as ss

    bamfh = prefix + sample_id + ".bam"
    bam = pysam.AlignmentFile(bamfh)
    
    ref = v.REF
    if isinstance(v.ALT, list):
        alt = v.ALT[0]
    else: alt = v.ALT

    strands = defaultdict(list)
    ref_mapq, alt_mapq = [], []

    for pileupcolumn in bam.pileup(v.CHROM, int(v.POS) - 1, int(v.POS) + 1):
        if pileupcolumn.pos != int(v.POS) - 1: continue
        for pileupread in pileupcolumn.pileups:
            # discard reads with gaps at site we're interested in
            if pileupread.is_del or pileupread.is_refskip: continue
            q_base = pileupread.alignment.query_sequence[pileupread.query_position]
            read_strand = '-' if pileupread.alignment.is_reverse else '+'
            if q_base.upper() == ref.upper():
                strands['ref'].append(read_strand) # strandedness of alignment
                ref_mapq.append(pileupread.alignment.mapping_quality) # MAPQ of alignment
            elif q_base.upper() == alt.upper():
                strands['alt'].append(read_strand)
                alt_mapq.append(pileupread.alignment.mapping_quality)
            else: continue
   
    ref_sum, alt_sum = len(strands['ref']), len(strands['alt'])
    total_depth = float(ref_sum + alt_sum)
    ab = alt_sum / total_depth if total_depth > 0 else 0.
    stat, p = ss.ranksums(ref_mapq, alt_mapq)
    mapq0_count = ref_mapq.count(0) + alt_mapq.count(0)

    aln_info = namedtuple('info', ['ab', 'mapq_stat', 'mapq_p', 'mapq0_count'])
    info = aln_info(ab, stat, p, mapq0_count)

    return info

def get_denovo(v, sample_dict, f1s, f2s,
        max_alts_in_parents=1,
        min_depth=5,
        max_mean_depth=400,
        min_allele_balance=0.1, 
        max_allele_balance=0.9,
        min_qual=1,
        min_p_ab=0.05, # Jonsson et al. (2017)
        max_het_outside_fam=10, 
        multigen=False,
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
    gts = v.gt_types
   
    # when we filter out potential cohort DNMs, we want
    # to look in everyone but the grandparents
    all_children = f1s + f2s

    if multigen: kids = f1s
    else: kids = f2s

    depths = None
    ref_depths = None
    # loop over alternate alleles
    for k in range(1, len(v.ALT) + 1):
        alt_depths = None
        # and then loop over the f1s. 
        for kid in kids:
            # Added check since sample '8477' apparently doesn't exist in our ceph data
            if simons_prefix and 'SSC' not in kid.sample_id: continue
            if kid.mom.sample_id == '8477' or kid.dad.sample_id == '8477':
                continue

            ki = sample_dict[kid.sample_id]
            mi, di = sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

            quals = v.gt_quals
            quals[quals < 0] == 0
            if quals[ki] < min_qual or quals[mi] < min_qual or quals[di] < min_qual: continue

            # basic sanity check for female f1s
            if v.CHROM == 'Y' and kid.sex == 'female': continue
            
            # check that kid is a HET (or HOM_ALT on hemizygous chroms in males)
            if kid.sex == 'male' and v.CHROM in ['X', 'Y']:
                if gts[ki] != 3: continue
            else:
                if gts[ki] != 1: continue
            
            # check that parents are homozygous
            # NOTE: what about parents that are HOM_ALT?
            if not (gts[di] == 0 and gts[mi] == 0) or (gts[di] == 3 and gts[mi] == 3): continue
            
            # check that likelihood of top genotype is 20 more than
            # likelihood of lowest genotype (from Jonsson et al. 2017)
            #gt_lh = sorted(v.format('PL')[ki])
            #if gt_lh[-1] < (gt_lh[-2] + 20): continue

            if depths is None:
                ref_depths = v.gt_ref_depths
                ref_depths[ref_depths < 0] = 0
            if alt_depths is None:
                alt_depths = v.gt_alt_depths
                alt_depths[alt_depths < 0] = 0
            
            k_alt, k_ref = alt_depths[ki], ref_depths[ki]
            m_alt, m_ref = alt_depths[mi], ref_depths[mi]
            d_alt, d_ref = alt_depths[di], ref_depths[di]
            
            # check total depth in mom, dad, and kid, and for lack of ALT support
            # in the parents
            if alt_depths[[mi, di]].sum() > max_alts_in_parents: continue
            if (m_ref + m_alt) < min_depth or (d_ref + d_alt) < min_depth: continue
            if (k_alt + k_ref) < min_depth: continue

            # strict allele balance filter in parents, from Jonsson et al. (2017)
            if m_alt / (m_alt + m_ref) > min_p_ab: continue
            if d_alt / (d_alt + d_ref) > min_p_ab: continue

            # if we're looking in a third generation, kids (F1s) have children
            children = []
            if multigen:
                children = [k for k in f2s if kid.sample_id in (k.dad.sample_id, k.mom.sample_id)]
            children_idxs = [sample_dict[k.sample_id] for k in children]
            # get the number of DNM genotypes in the cohort overall
            # NOTE: a male f1 with a "DNM" on the X chromosome is HOM_ALT, but
            # if a female in the cohort has the same "DNM," she's HET
            f1_dnm_outside_fam = 0 # real DNMs outside the fam in the F1
            f2_dnm_outside_fam = 0 # real DNMs outside the fam in the F2
            het_outside_fam = 0 # noise outside the fam
            gt_to_filter = 1
            if v.CHROM in ['X', 'Y'] and kid.sex == 'male': gt_to_filter = 3
            # get all VCF indexes where the person has the DNM genotype
            dnm_idxs = np.where(gts == gt_to_filter)[0]
            # exclude the index of the kid in question
            dnm_idxs_outside_kid = dnm_idxs[np.where(dnm_idxs != ki)[0]]
            # exclude the kid's children (if they have any)
            dnm_idxs_outside_fam = np.array([x for x in dnm_idxs_outside_kid if x not in children_idxs])
            try:
                # get GQ for all potential DNM
                dnm_gt_quals = quals[dnm_idxs_outside_fam]
                # get the potential DNMs above a GQ threshold
                dnm_above_qual = np.where(dnm_gt_quals > 20)[0]
                # get all F1 and F2 samples that have a potential DNM above the GQ threshold
                s_with_dnm_outside_fam = [s for s in all_children if sample_dict[s.sample_id] in dnm_above_qual] 
                for k in s_with_dnm_outside_fam:
                    if '8477' in (k.mom.sample_id, k.dad.sample_id): continue
                    m_idx, d_idx = sample_dict[k.mom.sample_id], sample_dict[k.dad.sample_id]
                    m_gt, d_gt = gts[m_idx], gts[d_idx]
                    # check if the sample is really a DNM (by genotype, that is)
                    if m_gt == d_gt and m_gt in (3, 0):
                        if k in f2s: f2_dnm_outside_fam += 1
                        elif k in f1s: f1_dnm_outside_fam += 1
                    # if not, it's just a spurious HET (or more likely, in the grandparents)
                    else: het_outside_fam += 1
            # if there aren't any DNM indices outside the family
            except IndexError: het_above_qual = 0 
            if f1_dnm_outside_fam + f2_dnm_outside_fam + het_outside_fam > max_het_outside_fam: continue
            
            # allele balance check in kid from VCF 
            ab = k_alt / (k_alt + k_ref)
            if not (min_allele_balance < ab < max_allele_balance): continue
           
            # check for inheritance of the DNM in a third generation
            f2_with_v = []
            p_dnm_in_f2s = 0.
            if multigen:
                f2_with_v = validate_in_f2(v, sample_dict, gts, kid, children)
                p_dnm_in_f2s = 1 - 0.5 ** (len(children) - len(f2_with_v))

            # impose lenient filters on individuals that inherit the DNM
            for f2 in f2_with_v:
                f2_idx = sample_dict[f2.sample_id]
                total_depth = alt_depths[f2_idx] + ref_depths[f2_idx]
                ab = alt_depths[f2_idx] / float(total_depth)
                if not (min_allele_balance < ab < max_allele_balance): f2_with_v.remove(f2)
                if total_depth < min_depth: f2_with_v.remove(f2) 
        
            # confirm that parents don't have excess of ALT reads in the BAM.
            # this is necessary because the GATK VCF often doesn't report all
            # ALT/REF support. also, check for an enrichment of bad alignments.
            ref_sb, alt_sb = 0., 0.
            mapq_stat, mapq_p = 0., 0.
            if check_bams:
                k_aln_info = query_alignments(v, kid.sample_id)
                d_aln_info = query_alignments(v, kid.mom.sample_id)
                m_aln_info = query_alignments(v, kid.mom.sample_id)
                if k_aln_info.mapq0_count > 5: continue
                if m_aln_info.ab > min_p_ab or d_aln_info.ab > min_p_ab: continue
                if not (min_allele_balance < k_aln_info.ab < max_allele_balance): continue

                mapq_stat, mapq_p = k_aln_info.mapq_stat, k_aln_info.mapq_p

            ret.extend(variant_info(v, kid, sample_dict, f2_with_v=f2_with_v, 
                                                        alt_i=k,
                                                        p_dnm_in_f2s=p_dnm_in_f2s,
                                                        ref_sb=ref_sb, alt_sb=alt_sb,
                                                        dnm_outside_fam=[str(f1_dnm_outside_fam), str(f2_dnm_outside_fam), str(het_outside_fam)],
                                                        mapq_stat=mapq_stat,
                                                        mapq_p=mapq_p))
            
    # shouldn't have multiple samples with same de novo.
    # TODO: make this a parameter. And/or handle sites where we get denovos from
    # multiple alts.
    if len(ret) == 1: return ret[0]

def variant_info(v, kid, sample_dict, f2_with_v=[], 
                    interval_type=None, alt_i=None, p_dnm_in_f2s=0.,
                    ref_sb=0., alt_sb=0., dnm_outside_fam=[],
                    mapq_stat=0., mapq_p=0.):

    ki, mi, di = sample_dict[kid.sample_id], sample_dict[kid.mom.sample_id], sample_dict[kid.dad.sample_id]

    f2_with_v_idx = [sample_dict[f2.sample_id] for f2 in f2_with_v]
    #depths = v.format('AD', int)
    #ref_depths = depths[:, 0]
    #all_alts = depths[:, 1:]

    ref_depths, alt_depths, quals = v.gt_ref_depths, v.gt_alt_depths, v.gt_quals

    k_ref, k_alt = ref_depths[ki], alt_depths[ki]

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
            ("dnm_outside_fam", ','.join(dnm_outside_fam)),
            ("num_hom_alt", v.num_hom_alt),
            ("allele_number", v.INFO.get("AN")),
            ("filter", v.FILTER or "PASS"),
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
            ("ref_strand_bias", ref_sb),
            ("alt_strand_bias", alt_sb),
            ("f2_with_variant", ','.join([str(x.sample_id) for x in f2_with_v])),
            ("p_dnm_in_f2s", p_dnm_in_f2s),
            ("mapq_ranksum_stat", mapq_stat),
            ("mapq_ranksum_p", mapq_p)
            ))

def denovo(v, sample_dict, f1s, f2s, max_alts_in_parents=2,
        min_depth=10, 
        max_mean_depth=400,
        min_allele_balance=0.1,
        max_allele_balance=0.9,
        min_p_ab=0.05,
        min_depth_percentile=3,
        min_qual=1,
        max_het_outside_fam=10,
        multigen=False,
        _use_cohort_filters=True,
        exclude=None,
        simons_prefix=False,
        check_bams=False,
        HET=1):
    if not variant_prefilter(v, 10): return None
    return get_denovo(v, sample_dict, f1s, f2s, max_alts_in_parents=max_alts_in_parents,
            min_depth=min_depth,
            max_mean_depth=max_mean_depth,
            min_allele_balance=min_allele_balance,
            max_allele_balance=max_allele_balance,
            min_qual=min_qual,
            max_het_outside_fam=max_het_outside_fam,
            multigen=multigen,
            min_p_ab=min_p_ab,
            exclude=exclude,
            _use_cohort_filters=_use_cohort_filters,
            simons_prefix=simons_prefix,
            check_bams=check_bams,
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

        d = denovo(v, sample_dict, f1s, f2s, min_depth=args.min_depth,
                min_allele_balance=args.min_allele_balance, 
                max_allele_balance=args.max_allele_balance, 
                max_alts_in_parents=args.max_alts_in_parents,
                min_p_ab=args.min_p_ab,
                min_qual=args.min_qual,
                max_het_outside_fam=args.max_het_outside_fam,
                multigen=args.multigen,
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
    p.add_argument("-min_allele_balance", help="min allele balance [0.1]", type=float, default=0.1)
    p.add_argument("-max_allele_balance", help="max allele balance [0.9]", type=float, default=0.9)
    p.add_argument("-max_alts_in_parents", help="max combined alt support in both parents [2]", type=int, default=2)
    p.add_argument("-min_qual", help="minimum GQ in kid [1]", type=int, default=1)
    p.add_argument("-min_p_ab", help="minimum allele balance allowed in parents [0.05]", type=float, default=0.05)
    p.add_argument("-max_het_outside_fam", help="maximum number of HET observations of the DNM outside the family [10]", type=float, default=10)
    p.add_argument("-simons_prefix", help="require all sample names to begin with the 'SSC' prefix", action="store_true")
    p.add_argument("-multigen", help="use a third generation of samples to validate DNMs", action="store_true")
    p.add_argument("-check_bams", help="check parental BAMs for evidence of DNM", action="store_true")
    p.add_argument("ped")
    p.add_argument("vcf")

    run(p.parse_args(argv))

if __name__ == "__main__":
    main(sys.argv[1:])
    if __package__ is None:
        from os import path
        sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
