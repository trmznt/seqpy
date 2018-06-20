
import math
#import pylab # this also import syshook to PyGTK stuff
import numpy
import colorsys

blv62 = {
    'A': ( -0.43, 4.02, 1.99 ),
    'C': ( 6.36, 5.21, 6.80),
    'L': ( 8.42, 3.01, -3.36)
}


hex_black = norm_black = ( 0, 0, 0 )
hex_white = ( 255, 255, 255)
norm_white = ( 1, 1, 1)

def bg_hex(r, g, b):
    if r + g + b > 255 * 1.25:
        return hex_black
    return hex_white

def bg_norm(r, g, b):
    if r + g + b > 1.25:
        return norm_black
    return norm_white

def hex2norm(r, g, b):
    return (1.0 * r / 255, 1.0 * g / 255, 1.0 * b / 255)

def norm2hex(r, g, b):
    return (int(r * 255), int(g * 255), int(b * 255))

class SeqColor(object):

    def __init__(self):
        self.hex = None
        self.norm = None

    def as_hex(self):
        return self.hex

    def as_norm(self):
        return self.norm

    def init_from_hex( self, colors ):
        colors_norm = self.convert_colors( colors, hex2norm )
        self.hex = self.prepare_contrast( colors, bg_hex)
        self.norm = self.prepare_contrast( colors_norm, bg_norm)


    def init_from_norm( self, colors ):
        values = []
        for c, r in colors.values():
            values += list( c )
            if r is not None:
                values += list( r )
        shift = 0 - min(values)
        highest = max(values) + shift

        colors_norm = self.convert_colors( colors,
                lambda r, g, b: ( (r+shift)/highest, (g+shift)/highest, (b+shift)/highest ) )
        colors_hex = self.convert_colors( colors_norm, norm2hex)
        self.norm = self.prepare_contrast( colors_norm, bg_norm )
        self.hex = self.prepare_contrast( colors_hex, bg_hex )



    def prepare_contrast(self, colors, func):
        d = {}
        for k in colors:
            c, r = colors[k]
            if r is None:
                r = func( *c )
            d[k] = c, r
        return d

    def convert_colors(self, colors, func):
        d = {}
        for k in colors:
            c, r = colors[k]
            c = func( *c )
            if r is not None:
                r = func( *r )
            d[k] = (c, r)
        return d


class HexSeqColor(SeqColor):

    colors = {}

    def __init__(self):
        super(SeqColor, self).__init__()
        self.init_from_hex( self.colors )

class NormSeqColor(SeqColor):

    colors = {}

    def __init__(self):
        super(SeqColor, self).__init__()
        self.init_from_norm( self.colors )
        

class EasyNA(HexSeqColor):

    colors = {
        'A': ( (182, 255, 159), None ),
        'T': ( (255, 159, 182), None ),
        'C': ( (159, 182, 255), None ),
        'G': ( (255, 220, 159), None ),
        'N': ( (0x33, 0x33, 0x33), None ),
        '-': ( (255, 255, 255), None ),
        '_unk_': ( (255, 255, 255), (128, 128, 128) )
    }


class AESSN3(NormSeqColor):

    colors = {
        'A': ( ( -0.99, -0.61,  0.00 ), None ),
        'R': ( (  0.28, -0.99, -0.22 ), None ),
        'N': ( (  0.77, -0.24,  0.59 ), None ),
        'D': ( (  0.74, -0.72, -0.35 ), None ),
        'C': ( (  0.34,  0.88,  0.35 ), None ),
        'Q': ( (  0.12, -0.99, -0.99 ), None ),
        'E': ( (  0.59, -0.55, -0.99 ), None ),
        'G': ( ( -0.79, -0.99,  0.10 ), None ),
        'H': ( (  0.08, -0.71,  0.68 ), None ),
        'I': ( ( -0.77,  0.67, -0.37 ), None ),
        'L': ( ( -0.92,  0.31, -0.99 ), None ),
        'K': ( ( -0.63,  0.25,  0.50 ), None ),
        'M': ( ( -0.80,  0.44, -0.71 ), None ),
        'F': ( (  0.87,  0.65, -0.53 ), None ),
        'P': ( ( -0.99, -0.99, -0.99 ), None ),
        'S': ( (  0.99,  0.40,  0.37 ), None ),
        'T': ( (  0.42,  0.21,  0.97 ), None ),
        'W': ( ( -0.13,  0.77, -0.90 ), None ),
        'Y': ( (  0.59,  0.33, -0.99 ), None ),
        'V': ( ( -0.99,  0.27, -0.52 ), None ),
        'X': ( (  0.99,  0.99,  0.99 ), None ),
        '-': ( (  0.25,  0.25,  0.25 ), None ),
        '_unk_': ( (0, 0, 0), None )
    }

def spectrum(size):
    step = 1.0/(size+1)
    return [ colorsys.hsv_to_rgb(x, 1, 0.75) for x in numpy.arange(0, 1.0, step) ]


