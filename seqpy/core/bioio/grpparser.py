

from seqpy import cout, cerr, cexit
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
    p.add_argument('--colourfile', default='')
    p.add_argument('--column', default='1,2')

    return p

# TODO: need to provide mechanism to select colour scheme

# 12 colours from ColorBrewer2
colour_list = [ '#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928',
                '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']

# 20 colours from Vega
colour_list = [ "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
                "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
                "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
                "#17becf", "#9edae5" ]

# 20 colours from IWantHue
colour_list = [ "#7ebbc7", "#6d38c0", "#72dc59", "#cc4dc4", "#cfd94a", "#6f68c9",
                "#649a40", "#c2477b", "#68d5a9", "#cd3f41", "#637dac", "#dc6931",
                "#4d285e", "#bf953c", "#cc9acd", "#536840", "#74372c", "#c9d19d",
                "#363638", "#c69085"]

class GroupParser(object):


    def __init__(self, args):


        self.groupfile = open(args.groupfile) if args.groupfile else None
        self.metafile = open(args.metafile) if args.metafile else None
        self.colourfile = open(args.colourfile) if args.colourfile else None
        self.delimiter = ( ',' if args.metafile[-3:].lower() == 'csv' else '\t' ) if args.metafile else None
        self.column = args.column
        self.group_info = {}        # { 'indv1': 'grp1', 'indv2': 'grp1', ...}
        self.groups = {}            # { 'grp1': [ 'indv1', 'indv2', ...], ...}
        self.group_keys = None      # [ 'grp1', 'grp1', 'grp2', 'grp1', ...]
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
            metadf = pandas.read_csv(self.metafile, sep=self.delimiter)
            sample_column, group_column = self.column.split(',')
            if sample_column.isdigit():
                sample_column = metadf.columns[ int(sample_column)-1 ]
            if group_column.isdigit():
                group_column = metadf.columns[ int(group_column)-1 ]

            cerr('[I: reading metafile for column: %s %s]'
                    % (sample_column, group_column))

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

        if not self.group_info:
            self.parse()

        groups = {}
        sample_idx = []
        group_keys = []
        for idx, code in enumerate(samples):
            grp_key = self.group_info[code]
            if grp_key in groups:
                groups[grp_key].append(idx)
            else:
                groups[grp_key] = [ idx ]
            sample_idx.append( idx )
            group_keys.append(grp_key)

        self.samples = samples
        self.sample_idx = set(sample_idx)
        self.groups = groups
        self.group_keys = group_keys

        if self.colourfile:
            # parse colour file
            self.colourfile.seek(0)
            next(self.colourfile)
            for line in self.colourfile:
                tokens = line.strip().split('\t')
                self.group_colours[tokens[0]] = tokens[1]

            # checking whether all groups has been assigned with colours
            for k in self.groups:
                if k not in self.group_colours:
                    cexit('E: group %s is not assigned' % k)

            cerr('[I: assigning manual colours to %d groups]' % (len(self.group_colours)))

        else:
            colour_wheel = cycle(colour_list)
            for k in sorted(self.groups.keys()):
                self.group_colours[k] = next(colour_wheel)

            if len(self.groups.keys()) > len(colour_list):
                cerr("W: warning, no of groups (%d) exceeds available colour list!" %
                    len(self.groups.keys()))

        return self.groups

    def colour_list(self, samples = None):
        """ return colour list based on samples """
        colours = []
        samples = samples or self.samples

        for s in samples:
            grp_key = self.group_info[s]
            colours.append( self.group_colours[grp_key] )

        return colours


    def group_colour_list(self, groups):
        """ return colour list based on groups """
        colours = []

        for g in groups:
            colours.append( self.group_colours[g] )

        return colours
