# tabutils.py

import pandas as pd
import pathlib

from enum import Enum

from seqpy import cerr

# tabular format for genomic data
#
# Sample-based Rows
# =================
#
# headers: SAMPLE, CHR:POS, ...
#
# each row contains a single sample
#
# Variant-based Column
#


def read_file(infile):
    match exts := pathlib.Path(infile).suffixes:
        case [*_, '.feather']:
            return pd.read_feather(infile)
        case [*_, '.txt'] | [*_, '.tsv'] | [*_, '.tsv', '.gz']:
            return pd.read_table(infile, sep='\t')
        case [*_, '.csv'] | [*_, '.csv', '.gz']:
            return pd.read_table(infile, sep=',')

    raise RuntimeError(f'Unknown extension type: {exts}')


def write_file(outfile, dataframe):
    match pathlib.Path(outfile).suffixes:
        case [*_, '.feather']:
            return dataframe.to_feather(outfile)
        case [*_, '.tsv'] | [*_, '.txt'] | [*_, '.tsv', '.gz'] | [*_, '.txt', '.gz']:
            return dataframe.to_csv(outfile, sep='\t', index=False)
        case [*_, '.csv'] | [*_, '.csv', '.gz']:
            return dataframe.to_csv(outfile, sep=',', index=False)
    raise RuntimeError('Unknown extension type')


def join_metafile(samples, metafile, percenttag=False):
    tokens = metafile.split(':')
    meta_df = read_file(tokens[0])
    if len(tokens) == 1:
        # we use all columns, and set the 1st columns as key
        columns = meta_df.columns
    else:
        columns = tokens[1].split(',')

        # if only 1 column name, assumed that it will be the indexing (ie sample code) column
        # and the all of the other columns will be joined
        if len(columns) == 1:
            columns = [columns[0]] + meta_df.columns

    joined_df = meta_df.meta.join_to_samples(samples, columns, percenttag)

    diff = None
    if len(joined_df) != len(samples):
        diff = set(samples) - set(joined_df['SAMPLE'])
        cerr('The following samples do not have metadata:')
        cerr(f'{diff}')

    return joined_df, diff


@pd.api.extensions.register_dataframe_accessor("meta")
class MetaAccessor(object):

    def __init__(self, df):
        self._df = df

    def join_to_samples(self, samples, columns, percenttag=False):

        actual_columns = []
        final_columns = []
        rename_dict = {}
        for c in columns:

            # if c contains '>', then the column will be renamed
            c, *c_rename = c.split('=', 1)

            # check if c exists
            if c not in self._df:
                raise ValueError(f'Column {c} does not exits')

            actual_columns.append(c)
            if any(c_rename):
                rename_dict[c] = c_rename[0]
                final_columns.append(c_rename[0])
            else:
                final_columns.append(c)

        # use pandas' join mechanism
        sample_df = pd.DataFrame(dict(SAMPLE=samples))
        joined_df = sample_df.join(self._df.loc[:, actual_columns].set_index(actual_columns[0]),
                                   on='SAMPLE', how='inner')

        # rename columns if needed
        if any(rename_dict):
            joined_df.rename(columns=rename_dict, inplace=True)

        # prepend percent signature if needed
        if percenttag:
            names = {}
            for c in final_columns[1:]:
                names[c] = f'%{c}'
            joined_df.rename(columns=names, inplace=True)

        return joined_df

    def join(self, right_df, on_column):

        # use pandas' join mechanism
        # import IPython; IPython.embed()
        return self._df.join(right_df.set_index(on_column), on=on_column, how='inner')

    def to_dict(self, key_column, value_column):
        return dict(zip(self._df[key_column], self._df[value_column]))


class FormatType(Enum):
    SAMPLE = 1
    VARIANT = 2


@pd.api.extensions.register_dataframe_accessor("geno")
class GenoAccessor(object):

    def __init__(self, df):
        self._type = None
        self._df = df
        self._sampledf = None
        self._metadf = None
        self._validate_process(df)

    def _validate_process(self, df):
        if df.columns[0] == 'SAMPLE':
            # sample-based rows
            self._type = FormatType.SAMPLE
            self._process_sample_format()

        elif df.columns[0] == 'CHROM' and df.columns[1] == 'POS':
            # position-based
            self._type = FormatType.VARIANT
            self._process_variant_format()

        else:
            raise AttributeError('This dataframe does not conform geno specs')

        return True

    def to_sampletype(self):

        # check that we start from position
        if self._type != 'variant':
            raise ValueError('Dataframe is already position-based')

    def to_positiontype(self):

        # check that we start from sample
        if self._type != 'sample':
            raise ValueError('Dataframe is already sample-based')

    def get_samplelist(self):
        if self._type == 'sample':
            return self._obj['SAMPLE']
        elif self._type == 'position':
            columns = self._obj.columns[2:]
            sample_indexes = ~columns.str.starswith('M|')
            return columns[sample_indexes]

    def get_positionlist(self):
        pass

    def get_positionids(self):
        pass

    def get_alleles(self):
        return self._alleledf

    def get_samples(self):
        return self._sampledf['SAMPLE']

    def _process_sample_format(self):
        # create new dataframe without metadata row but still with metadata columns
        metarows = self._df['SAMPLE'].str.startswith('%')
        allelecolumns = self._df.columns.str.contains(':')
        self._sampledf = self._df[~metarows]
        self._metarows = self._df[metarows]
        self._metacolumns = None
        self._alleledf = self._sampledf.loc[:, allelecolumns]

    def _process_variant_format(self):
        # create new dataframe without metadata columns
        raise NotImplementedError("the functionality is still being worked on")
        metacolumns = self._df.columns.str.startswith('%')
        self._metacolumns = self._df.loc[:, metacolumns]
        samplecolumns = self._df.columns[~metacolumns]

    def to_mhaps(self, mhap_specs):
        """ return a new DF containing microhaplotype strings """
        pass


@pd.api.extensions.register_dataframe_accessor("mhap")
class MicroHaplotypeAccessor(object):

    def __init__(self, df):
        self._df = df

    def to_numerics(self):
        """ return a new DF containing microhaplotype numerics """
        # for each mhap:
        #   perform value_count(), sort and convert to dict
        #   create new series based on the dict
        pass


def dataframe_from_variants(dataset, variants):

    from seqpy.core.sgk import sgutils

    position_ids = sgutils.get_position_ids(dataset)
    headers = ['SAMPLE'] + list(position_ids)

    columns = [dataset.sample_id] + variants

    df = pd.DataFrame(dict(zip(headers, columns)))

    return df


# EOF
