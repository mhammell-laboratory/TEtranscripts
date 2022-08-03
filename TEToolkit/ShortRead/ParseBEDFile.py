'''
Created on Oct 13, 2011

@author: Ying Jin
'''
import sys
import os
import struct
import logging
import struct
import gzip
from array import array
from time import time
from math import ceil,floor
from TEToolkit.Constants import BIN_SIZE,MAX_BIT, species_chrom_lengths
from TEToolkit.IO.FeatIO import FWTrackII
from TEToolkit.TEindex import *
#import rpy2.robjects as robjects


class Read :
    def __init__(self):
        self.chrom  = ""
        self.start  = -1
        self.end    = -1
        self.strand   = ""
        self.name = ""
        self.qual = 0

class GenericParser:
    """Generic Parser class.

    Inherit this to write your own parser.
    """
    def __init__ (self, filename):
        try :
            self.fhd = open(filename,'r')
        except:
            sys.exit(1)
        return

    def tsize(self):
        return

    def build_fwtrack (self):
        return

    def __fw_parse_line (self, thisline ):
        return

    def sniff (self):
        try:
            t = self.tsize()
        except:
            self.fhd.seek(0)
            return False
        else:
            if t<=10 or t>=10000:
                self.fhd.seek(0)
                return False
            else:
                self.fhd.seek(0)
                return t

class SAMFile :
    '''
      short reads in sam format
      saved by chromosome
   '''


    def __init__(self, srcfile,chroms):

        self.__srcfile = srcfile
        self.__fileList = dict()
        self.size = 0
        self.__binTags = []

class BEDFile(GenericParser) :
    '''
      short reads in bed format
      saved by chromosome
    '''


    def __init__(self, srcfile):
        '''
        Constructor
        '''
        self.__srcfile = srcfile
        self.__fileList = dict()
        self.size = 0
        self.__tsize = 0
        self.__binTags = []
        self.__buildAready = False

        self.fhd = open(srcfile, 'r')
        #self.__seprate_by_chrom(chroms)
        #self.__seprate_by_chrom()

    def sameFam(self,multi_reads,teIdx,seq_name=""):
        famlist = {}
        sel_reads = []
        #sel_idx = []

        for r in multi_reads :
            chr = r.chrom
            start = r.start-100
            end = r.start + 100
            famID = teIdx.getFamilyID(chr,start,end)
            if famID not in famlist :
                famlist[famID] = []
            famlist[famID].append(r)

        num_reads = len(multi_reads)
        max_fam1 = ""
        max_fam2 = ""
        max_fam1_cnt = 0
        max_fam2_cnt = 0
        for k in famlist:
            k_len = len(famlist[k])
            if max_fam1_cnt < k_len :
                max_fam1_cnt = k_len
                max_fam1 = k
            else :
                if max_fam2_cnt < k_len :
                    max_fam2_cnt = k_len
                    max_fam2_cnt = k
        if max_fam1_cnt > 2.5 * max_fam2_cnt :
            w = 1.0/len(famlist[max_fam1])
            return famlist[max_fam1],w
        else :
            sel_reads.extend(famlist[max_fam1])
            sel_reads.extend(famlist[max_fam2])
            w = 1.0/(len(famlist[max_fam1])+len(famlist[max_fam2]))
            return sel_reads, w

    def builtAready (self):
        return self.__buildAready

    def tsize(self):
        return self.__tsize

    def libsize(self):
        return self.size

    def build_fwtrack (self,temode):

        fwtrack = FWTrackII(filename=self.__srcfile)
        i = 0
        m = 0
        pre_seq_name = ""
        multi_reads = []
        seq_len_count = 0
        seq_len = 0
  #      cnt = 0
        strand = 0
        try:
            f = open(self.__srcfile,'r')
        except IOError:
            logging.error("open file %s error !\n" %(self.__srcfile))
            sys.exit(1)
        else:
            for line in f:
        # Go through bed file and assign each line to corresponding file.

                line = line.strip()
                items = line.split('\t')
                chrname = items[0]

                if seq_len_count < 1000 :
                    seq_len += int(items[2]) - int(items[1])
                    seq_len_count += 1

                i+=1
                start = 0
                if i == 1000000:
                        m += 1
                        logging.info(" %d" % (m*1000000))
                        i=0

                if items[5] == "+" :
                        strand = 0
                        start = int(items[1])
                else :
                        strand = 1
                        start = int(items[2])
                w = 1.0
                if len(items) > 6 : #there is weight assigned to each alignment
                    w = float(items[6])
                #    self.size += w
                    if temode == 'uniq' and w < 1.0 :
                       continue
                    fwtrack.add_loc(chrname,start,strand,w)

                else :
                    if pre_seq_name == "" :
                        pre_seq_name = items[3]

                    r = Read()
                    r.chrom = chrname
                    r.start = start
                    r.strand = strand
                    if pre_seq_name == items[3] :
                        multi_reads.append(r)
                    else :
                        if (temode =='uniq' and len(multi_reads) ==1) or (temode =='multi') :
                            w = 1.0/len(multi_reads)
                            for read in multi_reads :
                                fwtrack.add_loc(read.chrom,read.start,read.strand,w)
                            #self.size
                        multi_reads = []
                        pre_seq_name = items[3]
                        multi_reads.append(r)

              #  else:
               #     logging.warn("Unspecified chromosome name at %s line: %s. Skip!\n" %(self.__srcfile,line))

            if len(multi_reads) > 0 :
              if (temode =='uniq' and len(multi_reads) ==1) or (temode =='multi') :

                  w = float(1.0/len(multi_reads))
                  for k in range(len(multi_reads)) :
                     read = multi_reads[k]

                     fwtrack.add_loc(read.chrom,read.start,read.strand,w)

            f.close()
        if seq_len_count > 0 :
            fwtrack.setTsize(int(seq_len/seq_len_count))

        fwtrack.sort()
        self.__buildAready = True

        return fwtrack

    def build_fwtrack_v2 (self,teIdx):


        fwtrack = FWTrackII(filename=self.__srcfile)
        i = 0
        m = 0
        pre_seq_name = ""
        multi_reads = []
        seq_len_count = 0
        seq_len = 0
  #      cnt = 0
        strand = 0
        try:
            f = open(self.__srcfile,'r')
        except IOError:
            logging.error("open file %s error !\n" %(self.__srcfile))
            sys.exit(1)
        else:
            for line in f:
        # Go through bed file and assign each line to corresponding file.

                line = line.strip()
                items = line.split('\t')
                chrname = items[0]

                if seq_len_count < 1000 :
                    seq_len += int(items[2]) - int(items[1])
                    seq_len_count += 1

                i+=1
                start = 0

                if i == 1000000:
                        m += 1
                        logging.info(" %d" % (m*1000000))
                        i=0

                if items[5] == "+" :
                        strand = 0
                        start = int(items[1])
                        end = int(items[2])
                else :
                        strand = 1
                        start = int(items[2])

                w = 1.0
                if len(items) > 6 : #there is weight assigned to each alignment
                    w = float(items[6])
                #    self.size += w
                    fwtrack.add_loc(chrname,start,strand,w)

                else :
                    if pre_seq_name == "" :
                        pre_seq_name = items[3]

                    r = Read()
                    r.chrom = chrname
                    r.start = start

                    r.strand = strand
                    if pre_seq_name == items[3] :
                        multi_reads.append(r)
                    else :
                        (sel_reads,w) = self.sameFam(multi_reads,teIdx,pre_seq_name)
                       # w = 1.0/len(multi_reads)
                        for k in range(len(sel_reads)) :
                            read = sel_reads[k]


                            #w = weights[k]
                            fwtrack.add_loc(read.chrom,read.start,read.strand,w)
                            #self.size
                        multi_reads = []
                        pre_seq_name = items[3]
                        multi_reads.append(r)

              #  else:
               #     logging.warn("Unspecified chromosome name at %s line: %s. Skip!\n" %(self.__srcfile,line))

            if len(multi_reads) > 0 :

             # w = float(1.0/len(multi_reads))
              #for k in range(len(multi_reads)) :
               # read = multi_reads[k]
              (sel_reads,w) = self.sameFam(multi_reads,teIdx)
                       # w = 1.0/len(multi_reads)
              for k in range(len(sel_reads)) :
                            read = sel_reads[k]
                            #w = weights[k]
                            fwtrack.add_loc(read.chrom,read.start,read.strand,w)

            f.close()
        if seq_len_count > 0 :
            fwtrack.setTsize(int(seq_len/seq_len_count))

        fwtrack.sort()
        self.__buildAready = True

        return fwtrack


   # def __seprate_by_chrom (self,chroms):
    def __seprate_by_chrom (self):

        timestamp = time()
        fhead = "." + str(timestamp)

        #chroms = chrlen_tbl.keys()
        fhandels = {}
        seq_len = 0
        seq_len_count =1

        try:
            f = open(self.__srcfile,'r')
        except IOError:
            logging.error("open file %s error !\n" %(self.__srcfile))
            sys.exit(1)
        else:
            for line in f:
        # Go through bed file and assign each line to corresponding file.
                line = line.strip()
                items = line.split('\t')
                chrname = items[0]
                if seq_len_count < 1000 :
                    seq_len += int(items[2]) - int(items[1])
                    seq_len_count += 1
                if chrname not in fhandels :
                    chrfile = fhead + chrname + ".bed"
                    fh = open(chrfile,'a')
                    self.__fileList[chrname] = chrfile
                    fhandels[chrname] = fh

                    fh = fhandels[items[0]]

                    w = 1.0
                    if len(items) > 6 : #there is weight assigned to each alignment
                        w = float(items[6])
                    self.size += w
                    fh.write(items[0]+"\t"+items[1]+"\t"+items[2]+"\t"+items[3]+"\t"+items[4]+"\t"+items[5]+"\t"+str(w)+"\n")

                else :

                    fh = fhandels[items[0]]
                    #fh.write(line+"\n")
                    w = 1.0
                    if len(items) > 6 : #there is weight assigned to each alignment
                        w = float(items[6])
                    self.size += w
                    fh.write(items[0]+"\t"+items[1]+"\t"+items[2]+"\t"+items[3]+"\t"+items[4]+"\t"+items[5]+"\t"+str(w)+"\n")

              #  else:
               #     logging.warn("Unspecified chromosome name at %s line: %s. Skip!\n" %(self.__srcfile,line))

            f.close()
        self.__tsize = int(seq_len/seq_len_count)
        # Close all files.
        for fh in list(fhandels.values()):
            fh.close()

    def del_chrom_bed(self,chrom):
        f = self.__fileList[chrom]
        try:
            os.remove(f)

        except IOError :
            logging.error("cannot remove file %s !\n" %(f))
            sys.exit(1)

        #remove f

    def get_bin_rc (self,chrom,bin_size,fragsize,chrsize):

        reads = array('f',(0,) * int(ceil(1.0 *chrsize/bin_size)))
        sizeOfArray = len(reads)

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
                            loc = start + int(fragsize/2) -1
                        else :
                            loc = end - int(fragsize/2) -1
                        # Increment read count in bin.
                        pos = int(floor(1.0*loc/bin_size))

                        if(pos >= sizeOfArray) :
                            logging.error("%s  %d array index out of range! species may not match." %(items[0],pos))
                            sys.exit(0)
                        reads[pos] += w

                fh.close()
                return reads


class BAMFile(GenericParser) :
    '''
      short reads in BAM format
      saved by chromosome
      The bitwise flag is made like this:
    dec    meaning
    ---    -------
    1    paired read
    2    proper pair
    4    query unmapped
    8    mate unmapped
    16    strand of the query (1 -> reverse)
    32    strand of the mate
    64    first read in pair
    128    second read in pair
    256    alignment is not primary
    512    does not pass quality check
    1024    PCR or optical duplicate
    '''
    def __init__(self, srcfile):
        '''
        Constructor
        '''
        self.__srcfile = srcfile
        self.__fileList = dict()
#        self.size = 0
#        self.__tsize = 0

        self.__buildAready = False

        self.fhd = gzip.open(srcfile, 'r')
        try:
            self.fhd.read(10)
        except IOError:
        # not a gzipped file
            logging.error("it's not a BAM file %s !\n" %(srcfile))
            self.fhd.close()
            sys.exit(1)
        else:
            self.fhd.seek(0)

        #self.__seprate_by_chrom(chroms)
        #self.__seprate_by_chrom()

    def sameFam(self,multi_reads,teIdx):
        famlist = {}
        sel_reads = []
        #sel_idx = []

        for r in multi_reads :
            chr = r.chrom
            start = r.start - 100
            end = r.end + 100
            famID = teIdx.getFamilyID(chr,start,end)
            if famID not in famlist :
                famlist[famID] = []
            famlist[famID].append(r)

        num_reads = len(multi_reads)
        max_fam1 = ""
        max_fam2 = ""
        max_fam1_cnt = 0
        max_fam2_cnt = 0
        for k in famlist:
            k_len = len(famlist[k])
            if max_fam1_cnt < k_len :
                max_fam1_cnt = k_len
                max_fam1 = k
            else :
                if max_fam2_cnt < k_len :
                    max_fam2_cnt = k_len
                    max_fam2_cnt = k

        if max_fam1_cnt > 2.5 * max_fam2_cnt :
            w = 1.0/len(famlist[max_fam1])
            return famlist[max_fam1],w
        else :
            sel_reads.extend(famlist[max_fam1])
            sel_reads.extend(famlist[max_fam2])
            w = 1.0/(len(famlist[max_fam1])+len(famlist[max_fam2]))
            return sel_reads, w

    def builtAready (self):
        return self.__buildAready

    def build_fwtrack (self,temode):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII(filename=self.__srcfile)
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        references = []
        # move to pos 4, get the length of header
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]

        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            chrname = fread(nlength)[:-1]
            references.append(chrname)
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)

        i = 0
        m = 0
        multi_reads = []
        prev_seq = ""
        seq_len = 0
        seq_len_count = 0
        while 1:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:

                break
            (seq_name,chrid,fpos,flen,strand,qual) = self.__binary_parse(fread(entrylength))
            if seq_len_count < 1000 :
                seq_len += flen
                seq_len_count += 1


            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0

            if fpos >= 0:
                if seq_name == prev_seq : #multi reads
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = 1
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = 0
                    r.name = seq_name
                    r.weight = 0
                    r.qual = qual

                    multi_reads.append(r)
                else :
                    if prev_seq != "" :
                        if (temode == 'uniq' and len(multi_reads) == 1)  or temode == 'multi':
                            w = round(1.0/len(multi_reads),2)
                            for k in range(len(multi_reads)) :
                                rr = multi_reads[k]
                                fwtrack.add_loc(rr.chrom,rr.start,rr.strand,w)

                    multi_reads = []
                    prev_seq = seq_name
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = 1
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = 0
                    r.name = seq_name
                    r.qual = qual

                    multi_reads.append(r)


        if len(multi_reads) > 0 :
            if (temode == 'uniq' and len(multi_reads) == 1)  or temode == 'multi':
                w = round(1.0/len(multi_reads),2)
                for k in range(len(multi_reads)) :
                    rr = multi_reads[k]
                    fwtrack.add_loc(rr.chrom,rr.start,rr.strand,w)


        self.fhd.close()
        if seq_len_count > 0 :
            fwtrack.setTsize(int(seq_len/seq_len_count))

        self.__buildAready = True

        return fwtrack


    def build_fwtrack_v2 (self,teIdx):
        """Build FWTrackII from all lines, return a FWTrackII object.

        Note only the unique match for a tag is kept.
        """
        fwtrack = FWTrackII(filename=self.__srcfile)
        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        references = []
        # move to pos 4, get the length of header
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]

        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            chrname = fread(nlength)[:-1]
            references.append(chrname)
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)

        i = 0
        m = 0
        multi_reads = []
        prev_seq = ""
        seq_len = 0
        seq_len_count = 0
        while 1:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:

                break
            (seq_name,chrid,fpos,flen,strand,qual) = self.__binary_parse(fread(entrylength))
            if seq_len_count < 1000 :
                seq_len += flen
                seq_len_count += 1


            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0

            if fpos >= 0:
                if seq_name == prev_seq : #multi reads
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = 1
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = 0
                    r.name = seq_name
                    r.weight = 0
                    r.qual = qual

                    multi_reads.append(r)
                else :
                    if prev_seq != "" :
                        (sel_reads,w) = self.sameFam(multi_reads,teIdx)
                       # w = 1.0/len(multi_reads)
                        for k in range(len(sel_reads)) :
                            rr = sel_reads[k]
                            #w = weights[k]
                       # w = round(1.0/len(multi_reads),2)
                        #for k in range(len(multi_reads)) :
                         #   rr = multi_reads[k]
                            fwtrack.add_loc(rr.chrom,rr.start,rr.strand,w)

                    multi_reads = []
                    prev_seq = seq_name
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = 1
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = 0
                    r.name = seq_name
                    r.qual = qual

                    multi_reads.append(r)


        if len(multi_reads) > 0 :
            (sel_reads,w) = self.sameFam(multi_reads,teIdx)
                       # w = 1.0/len(multi_reads)
            for k in range(len(sel_reads)) :
                 rr = sel_reads[k]
                 #w = weights[k]
            #w = round(1.0/len(multi_reads),2)
            #for k in range(len(multi_reads)) :
             #   rr = multi_reads[k]
                 fwtrack.add_loc(rr.chrom,rr.start,rr.strand,w)


        self.fhd.close()
        if seq_len_count > 0 :
            fwtrack.setTsize(int(seq_len/seq_len_count))

        self.__buildAready = True

        return fwtrack


    def sniff(self):
        if self.fhd.read(3) == "BAM":
            return True
        else:
            return False

    #def __seprate_by_chrom (self,chroms):
    def __seprate_by_chrom (self):

        timestamp = time()
        fhead = "." + str(timestamp)

        #chroms = chrlen_tbl.keys()
        fhandels = {}

        # Create chromosome files for writing.
    #    for c in chroms:
    #        chrfile = fhead + c + ".bed"
    #        fh = open(chrfile,'a')
    #        self.__fileList[c] = chrfile
    #        fhandels[c] = fh

        fseek = self.fhd.seek
        fread = self.fhd.read
        ftell = self.fhd.tell
        references = []
        # move to pos 4, get the length of header
        fseek(4)
        header_len =  struct.unpack('<i', fread(4))[0]
        fseek(header_len + ftell())
        # get the number of chromosome
        nc = struct.unpack('<i', fread(4))[0]

        for x in range(nc):
            # read each chromosome name
            nlength = struct.unpack('<i', fread(4))[0]
            chrname = fread(nlength)[:-1]
            references.append(chrname)
            if chrname not in fhandels :
                chrfile = fhead + chrname + ".bed"
                fh = open(chrfile,'a')
                self.__fileList[chrname] = chrfile
                fhandels[chrname] = fh
            # jump over chromosome size, we don't need it
            fseek(ftell() + 4)

        i = 0
        m = 0
        multi_reads = []
        prev_seq = ""
        seq_len = 0
        seq_len_count = 0
        while 1:
            try:
                entrylength = struct.unpack('<i', fread(4))[0]
            except struct.error:

                break
            (seq_name,chrid,fpos,flen,strand,qual) = self.__binary_parse(fread(entrylength))
            if seq_len_count < 1000 :
                seq_len += flen
                seq_len_count += 1


            i+=1
            if i == 1000000:
                m += 1
                logging.info(" %d" % (m*1000000))
                i=0
               # print(fpos)
               # print(seq_name)
               # print(chrid)
               # print(flen)
            if fpos >= 0:
                if seq_name == prev_seq : #multi reads
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = '-'
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = '+'
                    r.name = seq_name
                    r.weight = 0
                    r.qual = qual

                    multi_reads.append(r)
                else :
                    if prev_seq != "" :
                        w = 1/len(multi_reads)
                        for k in range(len(multi_reads)) :
                            rr = multi_reads[k]
                            if rr.chrom in fhandels:
                                fh = fhandels[rr.chrom]
                                fh.write(rr.chrom+"\t"+str(rr.start)+"\t"+str(rr.end)+"\t"+rr.name+"\t"+str(rr.qual)+"\t"+rr.strand+"\t"+str(w)+"\n")
                                self.size += w
                            else :

                                logging.warn("Unspecified chromosome name at %s line: %s. Skip!\n" %(rr.chrom))
                    multi_reads = []
                    prev_seq = seq_name
                    r = Read()
                    r.chrom = references[chrid]
                    if strand == 1 :
                        r.strand = '-'
                        r.start = fpos - flen
                        r.end = fpos
                    else :
                        r.start = fpos
                        r.end = fpos + flen
                        r.strand = '+'
                    r.name = seq_name
                    r.qual = qual

                    multi_reads.append(r)


        if len(multi_reads) > 0 :
            w = 1/len(multi_reads)
            for k in range(len(multi_reads)) :
                rr = multi_reads[k]
                if rr.chrom in fhandels:
                    fh = fhandels[rr.chrom]
                    fh.write(rr.chrom+"\t"+str(rr.start)+"\t"+str(rr.end)+"\t"+rr.name+"\t"+str(rr.qual)+"\t"+rr.strand+"\t"+str(w)+"\n")
                    self.size += w

        self.fhd.close()
        self.__tsize = int(seq_len/seq_len_count)

        # Close all files.
        for fh in list(fhandels.values()):
            fh.close()

    def __binary_parse (self, data ):
        # we skip lot of the available information in data (i.e. tag name, quality etc etc)
        if not data: return (None,None,-1,-1,None,None)

        thisref = struct.unpack('<i', data[0:4])[0]
        thisstart = struct.unpack('<i', data[4:8])[0]
        (cigar, bwflag) = struct.unpack('<hh', data[12:16])
        (bin,mq_nl) = struct.unpack('<hh',data[8:12])
        name_len = mq_nl & 255

        seq_len = struct.unpack('<i',data[16:20])[0]
        seq_name = data[32:(32 + name_len)]
        seq_name = seq_name[:-1]
        qual = 0
        pairedEnd = False

        if bwflag & 4 or bwflag & 512 or bwflag & 1024:
            return (None,None, -1, -1,None,None)       #unmapped sequence or bad sequence
        if bwflag & 1:
            # paired read. We should only keep sequence if the mate is mapped
            # and if this is the left mate, all is within  the flag!
            if not bwflag & 2:
                return (None,None, -1, -1,None,None)   # not a proper pair
            if bwflag & 8:
                return (None,None, -1, -1,None,None)   # the mate is unmapped
            p1pos = thisstart
            p2pos = struct.unpack('<i', data[24:28])[0]
            if p1pos > p2pos:
                # this pair is the farthest one, skip it
                return (None,None, -1, -1,None,None)
        # In case of paired-end we have now skipped all possible "bad" pairs
        # in case of proper pair we have skipped the rightmost one... if the leftmost pair comes
        # we can treat it as a single read, so just check the strand and calculate its
        # start position... hope I'm right!
        if bwflag & 16:
            thisstrand = 1
            thisstart = thisstart + struct.unpack('<i', data[16:20])[0]    #reverse strand should be shifted len(query) bp
        else:
            thisstrand = 0

        return (seq_name,thisref, thisstart, seq_len,thisstrand,qual)


