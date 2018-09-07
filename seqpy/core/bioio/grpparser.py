

from seqpy import cout, cerr
from seqpy.cmds import arg_parser

from itertools import cycle

try:
    import pandas
except:
    cexit('ERR: rquire proper pandas instalation [pip3 install pandas]')

def init_argparser(p=None):

    if p is None:
        p = arg_parser('Group file parser')

    p.add_argument('--groupfile', default='')
    p.add_argument('--metafile', default='')
    p.add_argument('--column', default='1,2')

    return p

colour_list = [ '#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928',
                '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']

class GroupParser(object):


    def __init__(self, args):


        self.groupfile = open(args.groupfile) if args.groupfile else None
        self.metafile = open(args.metafile) if args.metafile else None
        self.delimiter = ( ',' if args.metafile[-3:].lower() == 'csv' else '\t' ) if args.metafile else None
        self.column = args.column
        self.groups = {}
        self.group_colours = {}



    def parse(self):

        # create a dictionary of groups <> sample_idx
        if self.groupfile:
            # this is a YAML/JSON file, open with YAML
            import yaml
            grouping = yaml.load( self.groupfile )
            groups = {}
            for g in grouping:
                for s in grouping[g]:
                    groups[s] = g
            self.group_info = groups

        elif self.metafile:
            # this is a tab/comma delimited file
            metadf = pandas.read_table(self.metafile, sep=self.delimiter)
            sample_column, group_column = self.column.split(',')
            if sample_column.isdigit():
                sample_column = metadf.columns[ int(sample_column)-1 ]
            if group_column.isdigit():
                group_column = metadf.columns[ int(group_column)-1 ]

            sampledf = metadf.loc[:, [sample_column, group_column] ]
            groups = {}
            for i in range(len(sampledf)):
                r = sampledf.loc[i]
                groups[r[0]] = r[1]

            self.group_info = groups

        else:
            cexit('E: need groupfile or metafile')

        return self.group_info


    def assign_groups(self, samples):

        groups = {}
        sample_idx = []
        for idx, code in enumerate(samples):
            grp_key = self.group_info[code]
            if grp_key in groups:
                groups[grp_key].append(idx)
            else:
                groups[grp_key] = [ idx ]
            sample_idx.append( idx )

        self.samples = samples
        self.sample_idx = set(sample_idx)
        self.groups = groups

        colour_wheel = cycle(colour_list)
        for k in sorted(self.groups.keys()):
            self.group_colours[k] = next(colour_wheel)
        return self.groups
