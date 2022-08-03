

"""Module for Feature IO classes.

Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@author:  Ying Jin
@contact: yjin@cshl.edu
"""

# ------------------------------------
# python modules
# ------------------------------------
import re
import sys
import logging
import struct
from array import array
from random import sample as random_sample
from operator import itemgetter
from math import sqrt,ceil,floor
import gzip
from time import time


from TEToolkit.Constants import *

# ------------------------------------
# Misc functions
# ------------------------------------


# ------------------------------------
# Classes
# ------------------------------------

class PeakIO:
    """IO for peak information.

    """
    def __init__ (self):
        self.peaks = {}

    def add (self, chromosome, start, end, summit=None,
             peak_height=None, number_tags=None,
             pvalue=None, fold_enrichment=None, fdr=None):
        """items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        if chromosome not in self.peaks:
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       peak_height,number_tags,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, pvalue_cut ):
        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut]
        self.peaks = new_peaks

    def filter_fc (self, fc_low, fc_up=None ):
        """Filter peaks in a given fc range.

        If fc_low and fc_up is assigned, the peaks with fc in [fc_low,fc_up)

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        chrs.sort()
        if fc_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low and p[7]<fc_up]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fc_low]
        self.peaks = new_peaks

    def total (self):
        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        x = 0
        for chrom in chrs:
            x += len(peaks[chrom])
        return x

    def ave_fold_enrichment (self):
        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        x = 0
        t = 0
        for chrom in chrs:
            x += len(peaks[chrom])
            for p in peaks[chrom]:
                t+=p[7]
        return t/x

    def max_fold_enrichment (self):
        """Return the maximum fc.

        """
        peaks = self.peaks
        chrs = list(peaks.keys())
        chrs.sort()
        x = 0
        for chrom in chrs:
            if peaks[chrom]:
                m = max([i[7] for i in peaks[chrom]])
                if m>x:
                    x=m
        return x


    def tobed (self):
        text = ""
        chrs = list(self.peaks.keys())
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (chrom,peak[0],peak[1])
        return text

    def towig (self):
        text = ""
        chrs = list(self.peaks.keys())
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text+= "%s\t%d\t%d\n" % (peak[0],peak[1])
        return text

    def init_from_dict (self, data):
        """Initialize the data from a dictionary. Improper assignment
        will damage the data structure.

        """
        self.peaks = {}
        chrs = list(data.keys())
        chrs.sort()
        for chrom in chrs:
            self.peaks[chrom]=[]
            a = self.peaks[chrom].append
            for i in data[chrom]:
                a(i)

    def merge_overlap ( self ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = list(peaks.keys())
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in range(len(peaks_chr)):
                if not prev_peak:
                    prev_peak = peaks_chr[i]
                    continue
                else:
                    if peaks_chr[i][0] <= prev_peak[1]:
                        s_new_peak = prev_peak[0]
                        e_new_peak = peaks_chr[i][1]
                        l_new_peak = e_new_peak-s_new_peak
                        if peaks_chr[i][4] > prev_peak[4]:
                            summit_new_peak = peaks_chr[i][3]
                            h_new_peak = peaks_chr[i][4]
                        else:
                            summit_new_peak = prev_peak[3]
                            h_new_peak = prev_peak[4]
                        prev_peak = (s_new_peak,e_new_peak,l_new_peak,summit_new_peak,h_new_peak,peaks_chr[i][5]+prev_peak[5])
                    else:
                        n_append(prev_peak)
                        prev_peak = peaks_chr[i]
            if prev_peak:
                n_append(prev_peak)
        del peaks
        self.peaks = new_peaks
        return True

    def overlap_with_other_peaks (self, peaks2, cover=0):
        """Peaks2 is a PeakIO object or dictionary with can be
        initialzed as a PeakIO. check __init__ for PeakIO for detail.

        return how many peaks are intersected by peaks2 by percentage
        coverage on peaks2(if 50%, cover = 0.5).
        """
        peaks1 = self.peaks
        if isinstance(peaks2,PeakIO):
            peaks2 = peaks2.peaks
        total_num = 0
        chrs1 = list(peaks1.keys())
        chrs2 = list(peaks2.keys())
        for k in chrs1:
            if not chrs2.count(k):
                continue
            rl1_k = iter(peaks1[k])
            rl2_k = iter(peaks2[k])
            tmp_n = False
            try:
                r1 = next(rl1_k)
                r2 = next(rl2_k)
                while (True):
                    if r2[0] < r1[1] and r1[0] < r2[1]:
                        a = sorted([r1[0],r1[1],r2[0],r2[1]])
                        if float(a[2]-a[1]+1)/r2[2] > cover:
                            if not tmp_n:
                                total_num+=1
                                tmp_n = True
                    if r1[1] < r2[1]:
                        r1 = next(rl1_k)
                        tmp_n = False
                    else:
                        r2 = next(rl2_k)
            except StopIteration:
                continue
        return total_num

class PosCnt:
    def __init__(self,pos,cnt):
        self.pos = pos
        self.cnt = cnt

    def __repr__(self):
        return repr((self.pos, self.cnt))


class FWTrackII:
    """Fixed Width Locations Track class II along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    def __init__ (self,filename = "", fw=0,anno=""):
        """fw is the fixed-width for all locations.

        """
        self.fw = fw
        self.__locations = {}           # locations
        self.__pluscounts = {}   # Ying read counts strand +
        self.__minuscounts = {}  # Ying read counts strand -
        self.__sorted = False
        self.total = 0                  # total tags
        self.__tagsize = 0
        self.__binTags = []
        self.annotation = anno   # need to be figured out
        self.__srcfile = filename

    def get_name(self):
        return self.__srcfile

    def clear_bins(self):
        del self.__binTags[:]

    def read_in_bins (self,species):

        size = 0.0
        idx_map = dict() #maps of chrom to index

        chrom_lengths = species_chrom_lengths[species]

       # for i in chrom_lengths.values() :
       #     size += i

        chrname = 'chr1'
        if species == 'dm3' :
            chrname = 'chr2L'

        size = chrom_lengths[chrname]

        for i in range(int(ceil(size/BIN_SIZE))):
            self.__binTags.append(0)

        pos = 0
        idx_map[chrname] = pos

        #read the bed file

                #    logging.warn("chrom %s does not exist in chrom_size file\n" %(c))
                #else :

        for poscnts in self.__pluscounts[chrname] :
            #start = poscnts.pos
            start = floor(poscnts)
            cnt = poscnts - start
            if start == poscnts :
                start -= 1
                cnt = 1
            offset = idx_map[chrname]+floor(start/BIN_SIZE)
            #self.__binTags[int(offset)] += poscnts.cnt
            self.__binTags[int(offset)] += cnt

        for poscnts in self.__minuscounts[chrname] :
            #start = poscnts.pos
            start = floor(poscnts)
            cnt = poscnts - start
            if start == poscnts :
                start -= 1
                cnt = 1
            offset = idx_map[chrname]+floor(start/BIN_SIZE)
            #self.__binTags[int(offset)] += poscnts.cnt
            self.__binTags[int(offset)] += cnt




    def get_bins(self,idx):

        cnts = array('f')
        pos = 0

        # idx is a byte array, each bit represent an element of self.__binTags

        for i in range(len(idx)):
            a = idx[i]
            pos = i * 8
            for j in range(8):
                b = a >> j
                if b & 1 ==1 :
                    cnts.append(1.0 * self.__binTags[pos+j])

        return cnts

    # return indices of bins with values fall in the fraction
    def get_bins_idx(self,fraction):

        idx = array('B')
        for i in range(MAX_BIT):
            idx.append(0x00)
       # qt  = robjects.r['quantile']
        #val = qt(robjects.FloatVector(self.__binTags),prob=fraction)

        for i in range(len(self.__binTags)):
            pos = int(floor(i/8))
            offset = i % 8
         #   if self.__binTags[i] <= val[0]  :
            if self.__binTags[i] > 0  :
                a = 1 << offset
                idx[pos] = idx[pos] | a

        return idx

    def get_all_bins_idx(self):

        idx = array('B')
        for i in range(MAX_BIT):
            idx.append(0x00)

        for i in range(len(self.__binTags)):
            pos = int(floor(i/8))
            offset = i % 8
            if self.__binTags[i] > 0 :
                a = 1 << offset
                idx[pos] = idx[pos] | a

        return idx

    def get_bin_rc_v2 (self,chrom,fragsize):

        #reads = array('f',(0,) * int(ceil(1.0 *chrsize/bin_size)))
        #sizeOfArray = len(reads)
        reads = []
        pos = []

        if chrom not in self.__fileList:
            logging.warn("No reads at chromosome %s. Skip!\n" %(chrom))
            return None
        else:
            fn = self.__fileList[chrom]
            try:
                fh = open(fn,'r')
            except IOError:
                logging.error("cannot open file of chromosome %s.\n" %(chrom))
                return None
            else:
                for line in fh:
                    line = line.strip()
                    items = line.split('\t')
                    if len(items) < 6 :
                        logging.warn("Format error at %s line: %s. Skip!\n" %(fn,line))
                    else:
                        start = int(items[1])
                        end = int(items[2])
                        strand = items[5]
                        w = 1
                        if len(items) > 6 : # weight for each tag
                            w = float(items[6])
                        # Shift read location by fragment size.
                        loc = 0
                        if strand == '+' :
                            loc = start + int(fragsize/2)
                        else :
                            loc = end - int(fragsize/2)
                        # Increment read count in bin.

                        #binIdx = int(floor(1.0*loc/bin_size))

                        #if(binIdx >= sizeOfArray) :
                        #    logging.error("%s  %d array index out of range! species may not match." %(items[0],pos))
                        #    sys.exit(0)
                        #reads[pos] += w

                        pos.append(loc)

                        reads.append(w)

                fh.close()
                return (reads,pos)


    def libsize(self):
        return self.total

    def add_loc (self, chromosome, fiveendpos, strand, cnt):
        """Add a location to the list according to the sequence name.

        chromosome -- mostly the chromosome name
        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """
        #if not self.__locations.has_key(chromosome):
        if chromosome not in self.__pluscounts :
            #self.__locations[chromosome] = [array(BYTE4,[]),array(BYTE4,[])] # for (+strand, -strand)
            self.__minuscounts[chromosome] = [] #[0]*100000 #[array(BYTE4,[])]
            self.__pluscounts[chromosome] = [] #[0]*100000 #[array(BYTE4,[])]
        #self.__locations[chromosome][strand].append(fiveendpos)
        #poscnt = PosCnt(fiveendpos,cnt)
        if strand == 0:  # Ying

            self.__pluscounts[chromosome].append(fiveendpos + cnt)
        else:
            self.__minuscounts[chromosome].append(fiveendpos + cnt)

        #self.total+=1

        self.total +=cnt # Ying total counts

    def get_locations_by_chr (self, chromosome):
        """Return a tuple of two lists of locations for certain chromosome.

        """

        if chromosome in self.__pluscounts:

            return (self.__pluscounts[chromosome],self.__minuscounts[chromosome])
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_locations_by_chr_v3 (self, chromosome):
        """Return a tuple of two lists of locations for certain chromosome.

        """

        if chromosome in self.__pluscounts:
            tagpos = []
            cnts = []
            pre_pos = -1
            cnt_sum = 0
            for i in range(len(self.__pluscounts[chromosome])) :
                poscnt = self.__pluscounts[chromosome][i]
                pos = int(floor(poscnt))
                cnt = poscnt - pos
                if pos == poscnt :
                    pos -= 1
                    cnt = 1
                if pre_pos == -1 :
                    pre_pos = pos
                    cnt_sum += cnt
                else :
                    if pre_pos == pos :
                        cnt_sum += cnt
                    else :
                        tagpos.append(pre_pos)
                        cnts.append(cnt_sum)
                        pre_pos = pos
                        cnt_sum = cnt
            if pre_pos != -1 :
                tagpos.append(pre_pos)
                cnts.append(cnt_sum)

            return (tagpos, cnts)
            #return (self.__pluscounts[chromosome],self.__pluscounts[chromosome])
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))

    def get_locations_by_chr_v2 (self, chromosome):
        """Return a tuple of two lists of locations for certain chromosome.

        """

        if chromosome in self.__pluscounts:
            tagpos = []
            cnts = []
            pre_pos = -1
           # cnt_sum = 0
            for i in range(len(self.__pluscounts[chromosome])) :
                poscnt = self.__pluscounts[chromosome][i]
                pos = int(floor(poscnt))
                cnt = poscnt - pos
                if pos == poscnt :
                    pos -= 1
                    cnt = 1

                tagpos.append(pos)
                cnts.append(cnt)

            return (tagpos, cnts)
            #return (self.__pluscounts[chromosome],self.__pluscounts[chromosome])
        else:
            raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))


    def clean(self,chr):
        self.__pluscounts[chr] = []

    def get_counts_by_chr (self, chromosome,strand=0): #Ying
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if strand ==0:
            if chromosome in self.__pluscounts:
                return self.__pluscounts[chromosome]
            else:
                raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))
        else:
            if chromosome in self.__minuscounts:
                return self.__minuscounts[chromosome]
            else:
                #raise Exception("No such chromosome name (%s) in TrackI object!\n" % (chromosome))
                return tmp;

    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        l = list(self.__pluscounts.keys())
        l.sort()
        return l

    def length (self):
        """Total sequenced length = total number of tags * width of tag
        """
        return self.total*self.fw

    def setTsize(self,tsize):
        self.__tagsize = tsize

    def tsize(self):
        return self.__tagsize

    def sort (self):
        """Naive sorting for locations.

        """
        for k in list(self.__pluscounts.keys()):

            (tmparrayplus,tmparrayminus) = self.get_locations_by_chr(k)
            self.__pluscounts[k] = sorted(tmparrayplus)
            self.__minuscounts[k] = sorted(tmparrayminus)
            #tmparray = self.get_counts_by_chr(k, 0) #Ying
            #self.__pluscounts[k] = sorted(tmparray,key=lambda poscnt: poscnt.pos)

            if len(tmparrayplus) < 1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
            #self.__locations[k][1] = sorted(tmparrayminus)
            #tmparray = []
            #tmparray = self.get_counts_by_chr(k, 1) #Ying
            #self.__minuscounts[k] = sorted(tmparray,key=lambda poscnt: poscnt.pos)
            if len(tmparrayminus) < 1:
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
        self.__sorted = True

    def filter_dup (self,maxnum):
        """Filter the duplicated reads.

        Run it right after you add all data into this object.
        """
        if not self.__sorted:
            self.sort()
        self.total = 0
        for k in self.__locations: # for each chromosome
            # + strand
            plus = self.__locations[k][0]
            if len(plus) <1:
                new_plus = []
            else:
                new_plus = array(BYTE4,[plus[0]])
                pappend = new_plus.append
                n = 1                # the number of tags in the current location
                current_loc = plus[0]
                for p in plus[1:]:
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            pappend(p)
                    else:
                        current_loc = p
                        pappend(p)
                        n = 1
                self.total +=  len(new_plus)

            # - strand
            minus = self.__locations[k][1]
            if len(minus) <1:
                new_minus = []
            else:
                new_minus = array(BYTE4,[minus[0]])
                mappend = new_minus.append
                n = 1                # the number of tags in the current location
                current_loc = minus[0]
                for p in minus[1:]:
                    if p == current_loc:
                        n += 1
                        if n <= maxnum:
                            mappend(p)
                    else:
                        current_loc = p
                        mappend(p)
                        n = 1
                self.total +=  len(new_minus)
            self.__locations[k]=[new_plus,new_minus]

    def merge_plus_minus_locations_naive (self):
        """Merge plus and minus strand locations

        """
        #for chrom in self.__locations.keys():
        for chrom in self.__pluscounts:
            #(plus_tags,minus_tags) = self.__locations[chrom]
            #self.__locations[chrom][0].extend(self.__locations[chrom][1])
            #self.__locations[chrom][0] = sorted(self.__locations[chrom][0])
            #self.__locations[chrom][1] = []
            #Ying
            self.__pluscounts[chrom].extend(self.__minuscounts[chrom])
            self.__pluscounts[chrom] = sorted(self.__pluscounts[chrom]) #,key=lambda x: x.pos)
            self.__minuscounts[chrom] = []



    def merge_plus_minus_locations (self):
        """Merge plus and minus strand locations.

        Tao: Amazingly, this function for merging two sorted lists is
        slower than merge_plus_minus_locations_naive which only
        concatenate the two lists then sort it again! I am so discouraged!
        """
        if not self.__sorted:
            self.sort()
        for chrom in self.__locations:
            (plus_tags,minus_tags) = self.__locations[chrom]
            new_plus_tags = array(BYTE4,[])
            ip = 0
            im = 0
            lenp = len(plus_tags)
            lenm = len(minus_tags)
            while ip < lenp and im < lenm:
                if plus_tags[ip] < minus_tags[im]:
                    new_plus_tags.append(plus_tags[ip])
                    ip += 1
                else:
                    new_plus_tags.append(minus_tags[im])
                    im += 1
            if im < lenm:
                # add rest of minus tags
                new_plus_tags.extend(minus_tags[im:])
            if ip < lenp:
                # add rest of plus tags
                new_plus_tags.extend(plus_tags[ip:])

            self.__locations[chrom] = [new_plus_tags,[]]
            self.total += len(new_plus_tags)


    def sample (self, percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        self.total = 0
        for key in self.__locations:
            num = int(len(self.__locations[key][0])*percent)
            self.__locations[key][0]=array(BYTE4,sorted(random_sample(self.__locations[key][0],num)))
            num = int(len(self.__locations[key][1])*percent)
            self.__locations[key][1]=array(BYTE4,sorted(random_sample(self.__locations[key][1],num)))
            self.total += len(self.__locations[key][0]) + len(self.__locations[key][1])

    def __str__ (self):
        return self.__to_wiggle()

    def __to_wiggle (self):
        """Use a lot of memory!

        """
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for k in self.__locations:
            if self.__locations[k][0]:
                t += "variableStep chrom=%s span=%d strand=0\n" % (k,self.fw)
                for i in self.__locations[k][0]:
                    t += "%d\t1\n" % i
            if self.__locations[k][1]:
                t += "variableStep chrom=%s span=%d strand=1\n" % (k,self.fw)
                for i in self.__locations[k][1]:
                    t += "%d\t1\n" % i
        return t

class FWTrackIII:
    """Fixed Width Locations Track class III along the whole genome
    (commonly with the same annotation type), which are stored in a
    dict.

    Locations are stored and organized by sequence names (chr names) in a
    dict. They can be sorted by calling self.sort() function.
    """
    def __init__ (self,fw=0,anno=""):
        """fw is the fixed-width for all locations.

        """
        self.fw = fw
        self.chr =""
        self.__locations = [array(BYTE4,[]),array(BYTE4,[])]         # locations
        self.__pluscounts = []   # Ying read counts strand +
        self.__minuscounts = []  # Ying read counts strand -
        self.__sorted = False
        self.total = 0                  # total tags
        self.annotation = anno   # need to be figured out

    def add_loc (self, fiveendpos, strand, cnt):
        """Add a location to the list according to the sequence name.


        fiveendpos -- 5' end pos, left for plus strand, right for neg strand
        strand     -- 0: plus, 1: minus
        """

        self.__locations[strand].append(fiveendpos)
        poscnt = PosCnt(fiveendpos,cnt)
        if strand == 0:  # Ying
            self.__pluscounts.append(poscnt)
        else:
            self.__minuscounts.append(poscnt)

        #self.total+=1
        self.total +=cnt # Ying total counts

    def get_locations(self):
        """Return a tuple of two lists of locations for certain chromosome.

        """

        return self.__locations


    def get_counts (self, strand): #Ying
        """Return a tuple of two lists of locations for certain chromosome.

        """
        if strand ==0:
            return self.__pluscounts[chromosome]
        else:
            return self.__minuscounts[chromosome]


    def get_chr_names (self):
        """Return all the chromosome names stored in this track object.
        """
        return self.chr

    def length (self):
        """Total sequenced length = total number of tags * width of tag
        """
        return self.total*self.fw

    def sort (self):
        """Naive sorting for locations.

        """
        for k in self.__locations:
            (tmparrayplus,tmparrayminus) = self.get_locations_by_chr(k)
            self.__locations[k][0] = sorted(tmparrayplus)
            tmparray = self.get_counts_by_chr(k, 0) #Ying
            self.__pluscounts[k] = sorted(tmparray,key=lambda poscnt: poscnt.pos)
            if len(tmparrayplus) < 1:
                logging.warning("NO records for chromosome %s, plus strand!" % (k))
            self.__locations[k][1] = sorted(tmparrayminus)
            tmparray = []
            tmparray = self.get_counts_by_chr(k, 1) #Ying
            self.__minuscounts[k] = sorted(tmparray,key=lambda poscnt: poscnt.pos)
            if len(tmparrayminus) < 1:
                logging.warning("NO records for chromosome %s, minus strand!" % (k))
        self.__sorted = True

    def merge_plus_minus_locations_naive (self):
        """Merge plus and minus strand locations

        """
        #(plus_tags,minus_tags) = self.__locations
        self.__locations[0].extend(self.__locations[1])
        self.__locations[0] = sorted(self.__locations[0])
        self.__locations[1] = []
            #Ying
        self.__pluscounts.extend(self.__minuscounts)
        self.__pluscounts.sorted(self.__pluscounts,key=lambda x: x.pos)
        self.__minuscounts = []

    def merge_plus_minus_locations (self):
        """Merge plus and minus strand locations.

        Tao: Amazingly, this function for merging two sorted lists is
        slower than merge_plus_minus_locations_naive which only
        concatenate the two lists then sort it again! I am so discouraged!
        """
        if not self.__sorted:
            self.sort()
        (plus_tags,minus_tags) = self.__locations
        new_plus_tags = array(BYTE4,[])
        ip = 0
        m = 0
        lenp = len(plus_tags)
        lenm = len(minus_tags)
        while ip < lenp and im < lenm:
           if plus_tags[ip] < minus_tags[im]:
               new_plus_tags.append(plus_tags[ip])
               ip += 1
           else:
               new_plus_tags.append(minus_tags[im])
               im += 1
        if im < lenm:
                # add rest of minus tags
          new_plus_tags.extend(minus_tags[im:])
        if ip < lenp:
               # add rest of plus tags
           new_plus_tags.extend(plus_tags[ip:])

        self.__locations = [new_plus_tags, []]
        self.total += len(new_plus_tags)


    def sample (self, percent):
        """Sample the tags for a given percentage.

        Warning: the current object is changed!
        """
        self.total = 0
        for key in self.__locations:
            num = int(len(self.__locations[key][0])*percent)
            self.__locations[key][0]=array(BYTE4,sorted(random_sample(self.__locations[key][0],num)))
            num = int(len(self.__locations[key][1])*percent)
            self.__locations[key][1]=array(BYTE4,sorted(random_sample(self.__locations[key][1],num)))
            self.total += len(self.__locations[key][0]) + len(self.__locations[key][1])

    def __str__ (self):
        return self.__to_wiggle()

    def __to_wiggle (self):
        """Use a lot of memory!

        """
        t = "track type=wiggle_0 name=\"tag list\" description=\"%s\"\n" % (self.annotation)
        for k in self.__locations:
            if self.__locations[k][0]:
                t += "variableStep chrom=%s span=%d strand=0\n" % (k,self.fw)
                for i in self.__locations[k][0]:
                    t += "%d\t1\n" % i
            if self.__locations[k][1]:
                t += "variableStep chrom=%s span=%d strand=1\n" % (k,self.fw)
                for i in self.__locations[k][1]:
                    t += "%d\t1\n" % i
        return t

