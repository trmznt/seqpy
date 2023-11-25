
import io

import pandas as pd
import numpy as np
import gzip

# This module is used to read files containing position in the genome.
#
# The files can be a BED-like format, which contain at least 3 corrdinate
# columns:
#   CHROM, START, END
# or a 2-column coordinate format:
#   CHROM, POS
#
# A file will be treated as BED-like file if:
# - its filename contains '.bed', or
# - its header contains CHROM START END
# otherwise, it will be treated as 2-column coordinate file.
#
# Arbitratry columns can follow after the coordinate columns
#


def init_argparser(p):
    p.add_argument('--posfile', default='')
    p.add_argument('--posfilefmt', default='')

    return p


def read_file(path):
    """ read file and discard lines starting with hash # """
    func = gzip.open if path.endswith('.gz') else open
    with func(path, 'rb') as fin:
        lines = []
        header = next(fin).strip()
        lines.append(header)
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith(b'#'):
                continue
            lines.append(line)
    buf = b'\n'.join(lines)
    return io.BytesIO(buf), header.decode('UTF-8')


def is_bedlike_format(filename, header):
    if '.bed' in filename.lower():
        return True

    tokens = header.split()
    if 'CHROM' in tokens and 'START' in tokens and 'END' in tokens:
        return True

    return False


def read_posfile(infile=None, args=None, sep=None, use_pyranges=False):
    """ read position file or bed-like file, and return plain pandas DataFrame
        or pyranges

        if infile contains '.bed', then treat as 2-coordinate bed file
    """

    if infile is None:
        if args:
            infile = args.posfile
        else:
            raise ValueError('read_posfile() need either a filename or args')

    buffer, header = read_file(infile)
    has_header = 0 if header.startswith('CHROM') else None
    df = pd.read_csv(buffer, sep=sep, header=has_header, engine='python')

    if is_bedlike_format(infile, header):
        # treat as 3-coordinate BED file
        if has_header is None:
            # add proper header to column 1, 2 and 3 as CHROM, START and END
            df.rename(columns={0: 'CHROM', 1: 'START', 2: 'END'}, inplace=True)

        # add _LENGTH column
        df['_LENGTH'] = df['END'] - df['START']

        # if all LENGTH is 1, create a POS column as well
        if df['_LENGTH'].eq(1).all():
            df['POS'] = df['END']

    else:
        if has_header is None:
            # add proper header to column 1 and 2 as CHROM, POS
            df.rename(columns={0: 'CHROM', 1: 'POS'}, inplace=True)
        # add _LENGTH column
        df['_LENGTH'] = 1
        df['START'] = df['POS'] - 1
        df['END'] = df['POS']
        # reorder columns to CHROM START END
        # code here

    if use_pyranges:
        try:
            import pyranges as pr
        except ModuleNotFoundError:
            raise RuntimeError('Please install pyranges: pip3 install pyranges')

        return pr.PyRanges(df)

    # by the end of the reading process, dataframe should have the following headers:
    # CHROM, START, END, _LENGTH and optional POS (if all _LENGTH == 1)

    return df


@pd.api.extensions.register_dataframe_accessor('pos')
class PositionAccessor(object):

    def __init__(self, df):
        self._validate(df)
        self._df = df

    @staticmethod
    def _validate(df):
        # check if at least have 2-column
        if 'CHROM' in df.columns and 'POS' in df.columns:
            if 'END' not in df.columns:
                df['START'] = df['POS'] - 1
                df['END'] = df['POS']
        else:
            raise ValueError('DataFrame does not have CHROM and POS')

    def point_tuples(self, as_list=False):
        # return a list of point position tuples [ (chrom, position), ...]
        if 'POS' in self._df.columns:
            if as_list:
                return (self._df['CHROM'], self._df['POS'])
            return zip(self._df['CHROM'], self._df['POS'])
        df_points = self._df[self._df['_LENGTH'] == 1]
        if as_list:
            return (df_points['CHROM'], df_points['END'])
        return zip(df_points['CHROM'], df_points['END'])

    def range_tuples(self):
        # return a list of range position tuples [ (chrom, start, end), ...]
        return zip(self._df['CHROM'], self._df['START'], self._df['END'])

    def get_indexes(self, source):
        """ return (indexes, targets)"""

        hash_table = {}

        # if source is a zarr dataset, get positions from it
        if (hasattr(source, 'contigs') and
                hasattr(source, 'variant_contig') and
                hasattr('variant_position')):
            snp_positions = list(
                    zip(np.array(source.contigs)[source.variant_contig.values],
                        source.variant_position.values))
            print(snp_positions)

        else:
            # assume that source is a list of (chrom, position, ...)
            snp_positions = source

        # create hash table to speed-up lookup
        for i, p in enumerate(snp_positions):
            hash_table[(p[0], p[1])] = i

        indexes = []
        targets = []

        for p in self.point_tuples():
            try:
                idx = hash_table[p]
                indexes.append(idx)
                targets.append(source[idx])
            except ValueError:
                indexes.append(-1)
                targets.append((p[0], p[1], None))

        return indexes, targets

    def sel_dataset(self, dataset):
        """ return a new dataset with variants based on positions """

        # convert position to (contig, positions)
        # currently only handles last position (END), ie. point position

        #import IPython; IPython.embed()

        contig_ids = dataset.contig_id.values.tolist()
        contigs = [contig_ids.index(c) for c in self._df['CHROM']]
        positions = list(zip(contigs, self._df['END']))

        return dataset.set_index(
            variants=('variant_contig', 'variant_position')
        ).sel(variants=positions).reset_index('variants').reset_coords(['variant_contig', 'variant_position'])

    def sel_tabular(self, tabframe):
        """ return a new tabular frame with variants based on positions """
        raise NotImplementedError()

    def to_bed(self, outpath, minimal=False):
        """ write to bedlike-file format """

        _df = self._df

        # make sure we have CHROM BEGIN END for the first 3 columns

        if not np.all(_df.columns[:3] == ['CHROM', 'START', 'END']):
            # rearrange columns
            for i, col in enumerate(['CHROM', 'START', 'END']):
                _df.insert(i, col, _df.pop(col))

        # remove unnecessary columns and create a copy to leave the original intact
        for col in ['_LENGTH', 'POS']:
            if col in _df.columns:
                _df = _df.drop(columns=col)

        if minimal:
            for col in _df.columns:
                if col in ['CHROM', 'START', 'END']:
                    continue
                _df = _df.drop(columns=col)

        # save to tab-delimited file
        _df.to_csv(outpath, sep='\t', index=False, header=False if minimal else True)

    def to_pos(self, outpath):
        """ write to position file format """

        _df = self._df

        # sanity check to ensure we don't have range position
        if np.any(_df['_LENGTH'] > 1):
            raise ValueError('to_pos() can only be used when all position have 1 bp length')

        # make sure we have CHROM POS for the first 2 columns
        if not np.all(_df.columns[:2] == ['CHROM', 'POS']):
            # rearange columns
            for i, col in enumerate(['CHROM', 'POS']):
                _df.insert(i, col, _df.pop(col))

        # remove unnecessary columns
        for col in ['_LENGTH', 'START', 'END']:
            if col not in _df.columns:
                continue
            _df = _df.drop(columns=col)

        # save to tab-delimited file
        _df.to_csv(outpath, sep='\t', index=False)

    def combine(self, df):
        """ merge between 2 position dataframe """

        # concatenate, but check for sanity first
        if ('POS' in self._df.columns) ^ ('POS' in df.columns):
            raise ValueError('Both must either have POS or not have POS column')

        _newdf = pd.concat([self._df, df])

        # if has POS column, sort by CHROM & POS and just drop duplicate straight away
        if 'POS' in _newdf.columns:
            return _newdf.drop_duplicates(subset=['CHROM', 'POS']).sort_values(by=['CHROM', 'POS'])

        # for range-based position, need to calculate overlap, possibly using pyranges
        raise NotImplementedError('The function has not been implemented yet')

    def intersection(self, df):
        """ create intersection between 2 position dataframe """

        # sanity check
        if ('POS' in self._df.columns) ^ ('POS' in df.columns):
            raise ValueError('Both must either have POS or not have POS column')

        if 'POS' in self._df.columns:
            keys = ['CHROM', 'POS']
            idx = pd.Index(self._df[keys]).intersection(df[keys])
            new_df = self._df.set_index(keys)
            return new_df[new_df.index.isin(idx)].reset_index().sort_values(by=keys)

        raise NotImplementedError("This functionality has not been implemented")

    def difference(self, df):
        """ create difference between this dataframe and df """

        # sanity check
        if ('POS' in self._df.columns) ^ ('POS' in df.columns):
            raise ValueError('Both must either have POS or not have POS column')

        if 'POS' in self._df.columns:
            keys = ['CHROM', 'POS']
            idx = pd.Index(self._df[keys]).difference(df[keys])
            new_df = self._df.set_index(keys)
            return new_df[new_df.index.isin(idx)].reset_index().sort_values(by=keys)

        raise NotImplementedError("This functionality has not been implemented")

    def symdiff(self, df):
        """ create symmetric difference between 2 position dataframe """

        # sanity check
        if ('POS' in self._df.columns) ^ ('POS' in df.columns):
            raise ValueError('Both must either have POS or not have POS column')

        if 'POS' in self._df.columns:
            keys = ['CHROM', 'POS']
            idx = pd.Index(self._df[keys]).symmetric_difference(df[keys])
            new_df = self._df.set_index(keys)
            return new_df[new_df.index.isin(idx)].reset_index().sort_values(by=keys)

    def to_mhcode(self):
        """ return a list of microhaplotype nomenclature code:
            CHROM:1ST_SNP:2ND,3RD,4TH
        """
        if 'MHAPIDX' not in self._df.columns:
            raise ValueError(
                'This method requires column named MHAPIDX for the microhaplotype index'
            )

        curr_chrom = curr_nomenclature = None
        curr_pos = curr_idx = -1
        nomenclatures = []
        ids = []

        for (idx, r) in self._df.iterrows():

            chrom, pos, idx = r['CHROM'], r['POS'], r['MHAPIDX']

            if curr_chrom is None:
                curr_chrom = chrom
                curr_pos = pos
                curr_idx = idx
                curr_nomenclature = f'{curr_chrom}:{curr_pos}:'
                continue

            if curr_idx != idx:
                nomenclatures.append(curr_nomenclature[:-1])
                ids.append(curr_idx)
                curr_chrom = chrom
                curr_pos = pos
                curr_idx = idx
                curr_nomenclature = f'{curr_chrom}:{curr_pos}:'
                continue

            curr_nomenclature = curr_nomenclature + f'{pos - curr_pos},'

        nomenclatures.append(curr_nomenclature[:-1])
        ids.append(curr_idx)

        return pd.DataFrame({'MHCODE': nomenclatures, 'ID': ids})

    # end of PositionAccessor()


def posframe_from_dataset(dataset, positions=None):
    """ create a position dataframe from a zarr dataset """

    # check type of positions
    match positions:
        case None:
            pass
        case list(tuple):
            # if a list of tuple
            # we need to translate chrom notation to numeric contig
            positions = list(zip(*positions))
            variant_contigs = [dataset.contigs.index(x) for x in positions[0]]
            positions = list(zip(variant_contigs, positions[1]))
            dataset = dataset.set_index(
                variants=('variant_contig', 'variant_position')
            ).sel(variants=positions)
        case list(str):
            # if a list of string (possibly CHROM:POS)
            positions = [p.split(':') for p in positions]
            positions = [(dataset.contigs.index(c), int(p)) for c, p in positions]
            dataset = dataset.set_index(
                variants=('variant_contig', 'variant_position')
            ).sel(variants=positions)
        case _:
            raise ValueError(f'Type {type(positions)} is not recognized for positions')

    columns = [
        ('CHROM', np.array(dataset.contigs)[dataset.variant_contig.values]),
        ('POS', dataset.variant_position.values),
        ('REF', dataset.variant_allele[:, 0]),
        ('ALT', dataset.variant_allele[:, 1]),
        ('_LENGTH', 1),
    ]

    if hasattr(dataset, 'variant_SNPEFF_GENE_NAME'):
        columns.append(('GENE', dataset.variant_SNPEFF_GENE_NAME.values))
    if hasattr(dataset, 'variant_SNPEFF_AMINO_ACID_CHANGE'):
        columns.append(('AACHANGE', dataset.variant_SNPEFF_AMINO_ACID_CHANGE))
    if hasattr(dataset, 'variant_SNPEFF_CODON_CHANGE'):
        columns.append(('CODON', dataset.variant_SNPEFF_CODON_CHANGE))

    return pd.DataFrame(dict(columns))


def posframe_from_tabular(dataframe):
    """ create a position dataframe from a tabular dataframe """
    raise NotImplementedError

# EOF
