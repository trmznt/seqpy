# tabutils.py

import pandas as pd


@pd.api.extensions.register_dataframe_accessor("geno")
class GenoAccessor:

    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    @staticmethod
    def _validate(obj):
        if obj.columns[0] == 'SAMPLE':
            # sample-based
            self._type = 'sample'

        elif obj.columns[0] == 'CHROM' and obj.columns[1] == 'POS':
            # position-based
            self._type = 'position'

        else:
            raise AttributeError('This dataframe does not conform geno specs')

        return True

    def to_sampletype(self):

        # check that we start from position
        if self._type != 'position':
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

    def get_allele_dataframe(self):
        pass


def dataframe_from_variants(dataset, variants):

    from seqpy.core.sgk import sgutils

    position_ids = sgutils.get_position_ids(dataset)
    headers = ['SAMPLE'] + list(position_ids)

    columns = [dataset.sample_id] + variants

    df = pd.DataFrame(dict(zip(headers, columns)))

    return df


# EOF
