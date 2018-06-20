

import cairo
import math


class SeqFigure(object):

    def __init__(self, mseq, positions, regions=None, linewidth=11, lineheight=14):
        self.mseq = mseq
        self.positions = positions
        self.regions = None
        self.linewidth = linewidth
        self.lineheight = lineheight


    def calc_max_label_width(self):
        pass

    def calc_max_pos_height(self):
        pass

    def calc_max_seq_width(self):
        pass

    def calc_dimension(self):
        pass


def draw_seqfig_at_surface( mseq, cairo_surface, colormaps):

    c = cairo.Context( cairo_surface )

    c.set_antialias( cairo.ANTIALIAS_NONE )



def calc_max_label_width( mseq,  ):

    max_width = 0

    for s in mseq:
        max_width = max( max_width, )