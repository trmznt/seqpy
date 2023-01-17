

# This library contains the necessary classes and methods to handle
# microhaplotype data. The nomenclature for the microhaplotype data
# is Chr:startpos:pos2,pos3,pos4 where pos2, pos3 etc are the offset
# of the consecutive SNPs from the first SNP position.
# eg: PvP01_01_v1:388122:41,75

import pandas as pd
from math import log, ceil
from numpy import base_repr
from seqpy import cerr
from seqpy.core.bioio import tabutils


def genotype_to_mhap(df, mhaplist, encode=False):
    """ Return a new pandas dataframe containing the microhaplotypes
        assembled from genotype dataframe, assuming that the dataframe
        contains single clonal sample or major genotype.
        Columns are [SAMPLE, CODE1, CODE2, ...]
    """

    sample_df = pd.DataFrame({'SAMPLE': df.geno.get_samples()})
    allele_df = df.geno.get_alleles()
    error_mhaps = []
    dfs = [sample_df]
    for mhcode in mhaplist:
        poslist = mhcode_to_columns(mhcode)

        try:
            geno_df = allele_df.loc[:, poslist]
        except KeyError:
            error_mhaps.append(mhcode)
            continue

        # sane checking
        if len(geno_df.columns) != len(poslist):
            raise ValueError(f'missing column for mhap: {mhcode}')

        # mark microhaplotypes that contain missing allele as missing mhap
        alleles = geno_df.iloc[:, 0].str.cat(geno_df.iloc[:, 1:], '')
        missing_idx = alleles.str.contains('X')
        alleles[missing_idx] = 'X'

        # mark microhaplotypes that contain heterozygous allele as missing mhap
        het_idx = alleles.str.contains('N')
        alleles[het_idx] = 'X'
        if het_idx.sum():
            cerr(f'[Warning: microhaplotypes found with heterozygous allele and '
                 f'marked as missing: {het_idx.sum()}]')

        dfs.append(
            pd.DataFrame({mhcode: alleles})
        )

    mhap_df = pd.concat(dfs, axis=1)

    return mhap_df, error_mhaps


def mhcode_to_columns(code):
    """ return a list of dataframe column id based from mhap code """
    chrom, pos, next_positions = code.split(':')
    pos = int(pos)
    next_positions = [int(p) for p in next_positions.split(',')]
    columns = [f'{chrom}:{pos}']
    for p in next_positions:
        columns.append(f'{chrom}:{pos + p}')

    return columns


# encoded haplotype uses base-5 integer

_to_int = {'-': '0', 'A': '1', 'C': '2', 'G': '3', 'T': '4'}
_to_na = {0: '-', 1: 'A', 2: 'C', 3: 'G', 4: 'T'}


def encode_haplotype(haplotype):
    alt_haplotype = ''.join(_to_int[a] for a in haplotype)
    return int(alt_haplotype, 5)


def decode_haplotype(numeric_hap, length):
    alt_haplotype = base_repr(numeric_hap, 5, length - ceil(log(numeric_hap, 5)))
    return ''.join(_to_na[a] for a in alt_haplotype)


# EOF
