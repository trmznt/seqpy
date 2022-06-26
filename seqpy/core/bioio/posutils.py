
import io

import pandas as pd
import numpy as np

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


def init_argparse(p):
    p.add_argument('--posfile', default='')
    p.add_argument('--posfilefmt', default='')

    return p


def read_file(path):
    """ read file and discard lines starting with hash # """
    with open(path) as fin:
        lines = []
        header = next(fin).strip()
        lines.append(header)
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                continue
            lines.append(line)
    buf = '\n'.join(lines)
    return io.StringIO(buf), header


def is_bedlike_format(filename, header):
    if '.bed' in filename.lower():
        return True

    tokens = header.split()
    if 'CHROM' in tokens and 'START' in tokens and 'END' in tokens:
        return True

    return False


def read_posfile(infile=None, args=None, use_pyranges=False):
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
    df = pd.read_csv(buffer, sep=None, header=has_header, engine='python')

    if is_bedlike_format(infile, header):
        # treat as 3-coordinate BED file
        if has_header is None:
            # add proper header to column 1, 2 and 3 as CHROM, START and END
            df.rename(columns={0: 'CHROM', 1: 'START', 2: 'END'}, inplace=True)
        # add _LENGTH column
        df['_LENGTH'] = df['END'] - df['START']

    else:
        if has_header is None:
            # add proper header to column 1 and 2 as CHROM, POS
            df.rename(columns={0: 'CHROM', 1: 'POS'})
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

    return df


@pd.api.extensions.register_dataframe_accessor('pos')
class PositionAccessor:

    def __init__(self, df):
        self._validate(df)
        self._df = df

    @staticmethod
    def _validate(df):
        # check if at least have 2-column
        pass

    def point_tuples(self):
        # return a list of point position tuples [ (chrom, position), ...]
        if 'POS' in self._df.columns:
            return zip(self._df['CHROM'], self._df['POS'])
        df_points = self._df[self._df['_LENGTH'] == 1]
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

        contigs = [dataset.contigs.index(c) for c in self._df['CHROM']]
        positions = list(zip(contigs, self._df['END']))

        return dataset.set_index(
            variants=('variant_contig', 'variant_position')
            ).sel(variants=positions)

    def sel_tabular(self, tabframe):
        """ return a new tabular frame with variants based on positions """
        raise NotImplementedError()

    def to_bed(self, outpath):
        """ write to bedlike-file format """

        _df = self._def

        # make sure we have CHROM BEGIN END for the first 3 columns

        if not np.all(_df.columns[:3] == ['CHROM', 'BEGIN', 'END']):
            # rearrange columns
            for i, col in enumerate(['CHROM', 'BEGIN', 'END']):
                _df.insert(i, col, _df.pop(col))

        # remove unnecessary columns and create a copy to leave the original intact
        for col in ['_LENGTH', 'POS']:
            if col in _df.columns:
                _df = _df.drop(columns=col)

        # save to tab-delimited file
        _df.to_csv(outpath, sep='\t', index=False)

    def to_pos(self, outpath):
        """ write to position file format """

        _df = self._df

        # sanity check to ensure we don't have range position
        if np.any(_df['_LENGTH' > 1]):
            raise ValueError('to_pos() can only be used when all position have 1 bp length')

        # make sure we have CHROM POS for the first 2 columns
        if not np.all(_df.columns[:2] == ['CHROM', 'POS']):
            # rearange columns
            for i, col in enumerate(['CHROM', 'POS']):
                _df.insert(i, col, _df.pop(col))

        # remove unnecessary columns
        for col in ['_LENGTH', 'BEGIN', 'END']:
            _df = _df.drop(columns=col)

        # save to tab-delimited file
        _df.to_csv(outpath, sep='\t', index=False)


def posframe_from_dataset(dataset):
    """ create a position dataframe from a zarr dataset """
    raise NotImplementedError()


def posframe_from_tabular(dataframe):
    """ create a position dataframe from a tabular dataframe """
    raise NotImplementedError

# EOF
