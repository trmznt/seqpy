
import numpy as np

from tqdm import tqdm


def _allele_for_barcode(datastore, variant_idx, allele_idxes, het_mask=None, failed_mask=None):

    alleles = datastore.variant_allele.sel(variants=variant_idx).values[allele_idxes]
    if het_mask is not None:
        alleles[het_mask] = 'N'
    if failed_mask is not None:
        alleles[failed_mask] = 'X'

    return alleles


def _allele_for_realmccoil(datastore, variant_idx, allele_idxes, het_mask=None, failed_mask=None):

    # for realmccoil, allele_idxes is the alleles
    # make sure we don't have idx > 1, otherwise convert to 1
    alleles = allele_idxes.astype(float)
    alleles[alleles > 1] == 1.0
    if het_mask is not None:
        alleles[het_mask] = 0.5
    if failed_mask is not None:
        alleles[failed_mask] = -1

    return alleles


def _allele_for_hmmibd(datastore, variant_idx, allele_idxes, het_mask=None, failed_mask=None):

    # for hmmibd, allele_idxes is the alleles
    alleles = allele_idxes
    if het_mask is not None:
        # hmmibd does not deal with het_mask, so just mark hets as missing
        alleles[het_mask] = -1
    if failed_mask is not None:
        alleles[failed_mask] = -1

    return alleles


def get_alleles(func, dataset, hetratio=0.67, mindepth=5, useGT=False):
    """ return a list of alleles per variant """

    variants = []
    ds = dataset
    for var_idx in tqdm(ds.variants):

        if useGT:

            genotypes = ds.call_genotype[var_idx]
            n_alt_alleles = genotypes.sum(axis=1)
            failed_mask = n_alt_alleles < 0
            het_mask = n_alt_alleles == 1

            # convert n_alt_alleles == 2 to 1 and treat as index
            n_alt_alleles[n_alt_alleles == 2] = 1

            alleles = func(ds, var_idx, n_alt_alleles, het_mask, failed_mask)

            variants.append(alleles)
            continue

        allele_depths = ds.call_AD[var_idx].values
        total_depths = allele_depths.sum(axis=1, where=(allele_depths > 0)) + 1e-3
        failed_mask = total_depths < mindepth

        # get allele based on highest depth
        allele_idx = np.argmax(allele_depths, axis=1)

        het_mask = None
        if hetratio > 0:
            # set heterozygosity based on depth ratio
            highest_depths = np.take_along_axis(
                allele_depths, np.expand_dims(allele_idx, axis=1), axis=1
                ).squeeze(axis=1)
            ratios = highest_depths / total_depths

            # use hetratio as cutoff for calling hets
            het_mask = ratios < hetratio

        alleles = func(ds, var_idx, allele_idx, het_mask, failed_mask)

        variants.append(alleles)

    return variants


def get_position_tuples(ds):
    return zip(np.array(ds.contigs)[ds.variant_contig.values], ds.variant_position.values)


def get_position_ids(ds):
    return [f'{c}:{p}' for c, p in get_position_tuples(ds)]

# EOF
