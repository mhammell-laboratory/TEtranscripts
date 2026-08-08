'''
Created on Oct 13, 2011

two normalization methods are supported.
sequence depth
bin correlation

@author: Ying Jin
'''
from array import array
import math
import struct
import zlib
#import numpy as np

import logging
import sys
from TEToolkit.Constants import *
from TEToolkit.ShortRead.ParseBEDFile import *


#def normalize(method,treatment,control,chrlen_tbl) :
def normalize(method,treatment,tinput,control,cinput,species,prj_name) :


    if method == 'sd' : #sequence depth
        return( seq_depth(treatment,tinput,control,cinput))
    elif method == 'bc' : #bin correlation
        #return( bin_corr(treatment,control,chrlen_tbl))
        return( bin_corr(treatment,tinput,control,cinput,species,prj_name))


def seq_depth(list1,list1input,list2,list2input):

    max_size = 0.0

    for i in range(len(list1)):
        tsmp = list1[i]
        if tsmp.libsize() > max_size :
            max_size = tsmp.libsize()


    for i in range(len(list2)):

        if list2[i].libsize() > max_size :
            max_size = list2[i].libsize()

    if len(list1input) > 0 and  max_size < list1input[0].libsize() :
        max_size = list1input[0].libsize()

    if len(list2input) > 0 and max_size < list2input[0].libsize() :
        max_size = list2input[0].libsize()

    sf1 = []
    sf2 = []

    for i in range(len(list1)):
        sf1.append(1.0 * max_size/list1[i].libsize())

    if len(list1input) > 0 :
        sf1.append(1.0*max_size/list1input[0].libsize())
    for i in range(len(list2)):
        sf2.append(1.0 * max_size/list2[i].libsize())
    if len(list2input) > 0 :
        sf2.append(1.0*max_size/list2input[0].libsize())

    return (sf1,sf2)

def _linear_scale_factor(reference, sample):
    """Return the zero-intercept OLS slope used by the former R helper."""

    denominator = sum(float(value) * float(value) for value in sample)
    if denominator == 0:
        return 1.0
    numerator = sum(float(left) * float(right) for left, right in zip(reference, sample))
    return numerator / denominator


def _png_chunk(chunk_type, data):
    payload = chunk_type + data
    return struct.pack(">I", len(data)) + payload + struct.pack(">I", zlib.crc32(payload) & 0xffffffff)


def _write_scatter_png(filename, reference, sample, width=500, height=500):
    """Write the historical log/log normalization scatter plot without R."""

    pixels = bytearray([255]) * (width * height * 3)

    def set_pixel(x_value, y_value, color):
        if 0 <= x_value < width and 0 <= y_value < height:
            offset = (y_value * width + x_value) * 3
            pixels[offset:offset + 3] = bytearray(color)

    left, right, top, bottom = 48, width - 16, 16, height - 42
    for x_value in range(left, right + 1):
        set_pixel(x_value, bottom, (50, 50, 50))
    for y_value in range(top, bottom + 1):
        set_pixel(left, y_value, (50, 50, 50))

    points = [
        (math.log10(float(x_value)), math.log10(float(y_value)))
        for y_value, x_value in zip(reference, sample)
        if x_value > 0 and y_value > 0
    ]
    if points:
        x_values = [point[0] for point in points]
        y_values = [point[1] for point in points]
        x_min, x_max = min(x_values), max(x_values)
        y_min, y_max = min(y_values), max(y_values)
        if x_min == x_max:
            x_min, x_max = x_min - 0.5, x_max + 0.5
        if y_min == y_max:
            y_min, y_max = y_min - 0.5, y_max + 0.5
        for x_log, y_log in points:
            x_pixel = int(round(left + (x_log - x_min) * (right - left) / (x_max - x_min)))
            y_pixel = int(round(bottom - (y_log - y_min) * (bottom - top) / (y_max - y_min)))
            for x_offset in (-1, 0, 1):
                for y_offset in (-1, 0, 1):
                    set_pixel(x_pixel + x_offset, y_pixel + y_offset, (38, 91, 160))

    scanlines = bytearray()
    stride = width * 3
    for row in range(height):
        scanlines.append(0)
        scanlines.extend(pixels[row * stride:(row + 1) * stride])

    signature = b"\x89PNG\r\n\x1a\n"
    header = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    with open(filename, "wb") as handle:
        handle.write(signature)
        handle.write(_png_chunk(b"IHDR", header))
        handle.write(_png_chunk(b"IDAT", zlib.compress(bytes(scanlines), 9)))
        handle.write(_png_chunk(b"IEND", b""))

#def bin_corr(list1,list2,chrlen_tbl):
def bin_corr(list1,list1input,list2,list2input,species,prj_name):

    tisize = len(list1input)
    cisize = len(list2input)
    tsize = len(list1)
    csize = len(list2)
    fnames = []

    #fnames2 = []
    sf1 = []
    sf2 = []


    for i in range(tsize):
        fnames.append(list1[i].get_name())
        list1[i].read_in_bins(species) # read  tags of treatment samples into bins

    if tisize > 0 :
        fnames.append(list1input[0].get_name())
        list1input[0].read_in_bins(species)



    for i in range(csize):
        fnames.append(list2[i].get_name())
       # list2[i].read_in_bins(chrlen_tbl) # read tags of control samples into bins
        list2[i].read_in_bins(species)

    if cisize > 0 :
        fnames.append(list2input[0].get_name())
        list2input[0].read_in_bins(species)

    # merge bins to make a matrix, rows are coordinates and cols are smps
    # randomly select the first sample as a temporary reference
    (tot_reads,reads,sel_idx) = join_bins(list1,list1input,list2,list2input,fnames,0)

    for i in range(csize):
        list2[i].clear_bins()

    if cisize > 0 :
        list2input[0].clear_bins()
    if tisize > 0 :
        list1input[0].clear_bins()

    for i in range(tsize):
        list1[i].clear_bins()

    # Calculate the same zero-intercept regression slope as ``lm(ref ~
    # sample - 1)`` and emit the same project-prefixed PNG artifacts natively.
    for i in range(tsize+tisize):
        if i == 0:
            msf = 1.0
        else:
            msf = _linear_scale_factor(reads[0], reads[i])
            outfname = prj_name+"_"+fnames[0]+"vs"+fnames[i]+".png"
            _write_scatter_png(outfname, reads[0], reads[i])
        sf1.append(msf)

    for i in range(tsize+tisize,csize + cisize + tsize+tisize):
        outfname = prj_name+"_"+fnames[0]+"vs"+fnames[i]+".png"
        msf = _linear_scale_factor(reads[0], reads[i])
        _write_scatter_png(outfname, reads[0], reads[i])
        sf2.append(msf)


    #plot scatter plot
    #sf = plot_sf(tot_reads,reads,sel_idx,fnames,min_sf_idx)

    #return (sf1[0:tsize],sf2[tsize:(tsize+csize)])
    return (sf1,sf2)


def join_bins(list1,list1input,list2,list2input,fnames,ref_idx):

    # only part of the background bins are used for computing normalization factors
    fraction = 1.0 #0.75
    idx_merged = array('B') #use a bit array to reduce memory usage
    reads = [] #np.array(dtype='double')
    tot_reads = []
    sel_idx = []
    #init a bit array
    for i in range(MAX_BIT):
        idx_merged.append(0x00)

    # union of indices that have reads >0
    for i in range(len(list1)):
        idx = list1[i].get_all_bins_idx()
        for j in range(len(idx)):
            idx_merged[j] = idx_merged[j] | idx[j]

    if len(list1input) > 0 :
        idx = list1input[0].get_all_bins_idx()
        for j in range(len(idx)):
            idx_merged[j] = idx_merged[j] | idx[j]


    for i in range(len(list2)):
        idx = list2[i].get_all_bins_idx()
        for j in range(len(idx)):
            idx_merged[j] = idx_merged[j] | idx[j]

    if len(list2input) > 0 :
        idx = list2input[0].get_all_bins_idx()
        for j in range(len(idx)):
            idx_merged[j] = idx_merged[j] | idx[j]

    #get total reads
#    for i in range(len(list1)):
#        r = list1[i].get_bins(idx_merged)
#        tot_reads.append(r)

#    if len(list1input) > 0 :
#        r = list1input[i].get_bins(idx_merged)
#        tot_reads.append(r)

#    for i in range(len(list2)):
#        r = list2[i].get_bins(idx_merged)
#        tot_reads.append(r)

#    if len(list2input) > 0 :
#        r = list2input[i].get_bins(idx_merged)
#        tot_reads.append(r)

    idx_tot = []
    idx_tot.extend(idx_merged)
    #intersect bins that are in the selected fraction
#    for i in range(len(list1)):
#        idx = list1[i].get_bins_idx(fraction)
#        for j in range(len(idx)):
#            idx_merged[j] = idx_merged[j] & idx[j]


#    for i in range(len(list2)):
#        idx = list2[i].get_bins_idx(fraction)
#        for j in range(len(idx)):
#            idx_merged[j] = idx_merged[j] & idx[j]

#    for i in range(len(idx_tot)):
#        a = idx_tot[i]
#        b = idx_merged[i]
#        if a | b != 0 :
#            if a & 1 == 1 and b & 1 == 1:
#                sel_idx.append(1)
#            if a & 1 == 1 and b & 1 == 0:
#                sel_idx.append(0)
#            for j in range(7):
#                a = a >> 1
#                b = b >> 1
#                if a & 1 == 1 and b & 1 == 1:
#                    sel_idx.append(1)
#                if a & 1 == 1 and b & 1 == 0:
#                    sel_idx.append(0)


    # retrieve reads in the bins (intersections of all samples)
    for i in range(len(list1)):
        r = list1[i].get_bins(idx_merged)
        reads.append(r)


    if len(list1input) > 0 :
        r = list1input[0].get_bins(idx_merged)
        reads.append(r)

    for i in range(len(list2)):
        r = list2[i].get_bins(idx_merged)
        reads.append(r)

    if len(list2input) >0 :
        r = list2input[0].get_bins(idx_merged)
        reads.append(r)

    if len(reads) == 0 :
        logging.error("normalization failed !")
        sys.exit(0)

    return (tot_reads,reads,sel_idx)
