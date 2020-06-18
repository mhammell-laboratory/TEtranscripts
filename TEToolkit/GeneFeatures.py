"""Module Description

    Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


    This code is free software; you can redistribute it and/or modify it
    under the terms of the Artistic License (see the file COPYING included
    with the distribution).

    @author:  Ying Jin
    @contact: yjin@cshl.edu
    """
import sys, time,re
import logging
import gzip
from math import ceil,floor
import collections

from TEToolkit.IntervalTree import *

#Taken from HTSeq
class GFF_Reader( ):

   """Parse a GFF file

   Pass the constructor either a file name or an iterator of lines of a
   GFF files. If a file name is specified, it may refer to a gzip compressed
   file.

   Yields tuple of (gene_id,chrom,strand,start position,end position,type)

   """

   def __init__( self, filename, id_attribute):
      self.line_no = None
      self.filename = filename
      self.id_attribute = id_attribute
      self._re_attr_main = re.compile( "\s*([^\s\=]+)[\s=]+(.*)" )

   def __iter__( self ):
      self.line_no = 0
      if self.filename.lower().endswith( ( ".gz" , ".gzip" ) ):
            lines = gzip.open( self.filename )
      else:
            lines = open( self.filename )

      for line in lines:
          self.line_no += 1
          if line == "\n" or line.startswith('#'):
              continue
          ( seqname, source, feature, start, end, score, strand, frame, attributeStr ) = line.split("\t")
          id = self.__parse_GFF_attr_string(attributeStr,self.id_attribute)

          yield (id, seqname, strand, int(start), int(end), feature)

      lines.close()
      self.line_no = None

   def __parse_GFF_attr_string(self,attributeStr,id_interested) :


       for pairs in attributeStr.split(';') :
           if pairs.count('"') not in [0,2] :
               raise ValueError("The attribute string seems to contain mismatched quotes.")
           nv = self._re_attr_main.match(pairs)
           if not nv :
               raise ValueError("Failure parsing GFF attribute line.")
           val = nv.group(2)
           name = nv.group(1)
           if name == id_interested :
               return val
       return None

   def get_line_number_string( self ):
      if self.line_no is None:
            return "file %s closed" % self.filename

      else:
         return "line %d of file %s" % ( self.line_no, self.filename )


class GeneFeatures:
    """index of Gene annotations.
        """
    def __init__ (self,GTFfilename,stranded,feature_type,id_attribute):

        self.featureIdxs_plus = {}
        self.featureIdxs_minus = {}
        self.featureIdxs_nostrand = {}
        self.features = []

        self.read_features(GTFfilename,stranded,feature_type,id_attribute)


    # Reading & processing annotation files
    def read_features(self,gff_filename, stranded, feature_type, id_attribute) :

        #dict of dicts since the builtin type doesn't support it for some reason
        temp_plus = collections.defaultdict(dict)
        temp_minus = collections.defaultdict(dict)
        temp_nostrand = collections.defaultdict(dict)

        # read count of features in GTF file
        gff = GFF_Reader(gff_filename,id_attribute)  # (id, seqname, strand, int(start), int(end), feature)
        i = 0
        counts = 0
        try:
            for f in gff:
                if f[0] is None :
                    continue
                if f[5] == feature_type:
                    counts += 1
                    if stranded != "no" and f[2] == "." :
                        sys.stderr.write("Feature %s does not have strand information." % (f[0]))
                    try:
                        if f[2] == "."  :
                            temp_nostrand[f[1]][f[0]].append((f[3],f[4]))
                    except:
                        temp_nostrand[f[1]][f[0]] = [(f[3],f[4])]

                    try:
                        if f[2] == "+"  :
                            temp_plus[f[1]][f[0]].append((f[3],f[4]))
                    except:
                        temp_plus[f[1]][f[0]] = [(f[3],f[4])]

                    try:
                        if f[2] == "-" :
                            temp_minus[f[1]][f[0]].append((f[3],f[4]))
                    except KeyError:
                        temp_minus[f[1]][f[0]] = [(f[3],f[4])]

                    #save gene id
                    if f[0] not in self.features :
                        self.features.append(f[0])

                    i += 1
                    if i % 100000 == 0 :
                        sys.stderr.write("%d GTF lines processed.\n" % i)
        except:
            sys.stderr.write("Error occured in %s.\n" % gff.get_line_number_string())
            raise

        if counts == 0 :
            sys.stderr.write("Warning: No features of type '%s' found in gene GTF file.\n" % feature_type)

        #build interval trees
        for each_chrom in temp_plus:
            inputlist = []
            for each_gene in temp_plus[each_chrom]:
                    for (start,end) in temp_plus[each_chrom][each_gene]:
                        inputlist.append(Interval(each_gene,start,end))
            self.featureIdxs_plus[each_chrom] = IntervalTree(inputlist)


        for each_chrom in temp_minus:
            inputlist = []
            for each_gene in temp_minus[each_chrom]:
                for (start,end) in temp_minus[each_chrom][each_gene]:
                    inputlist.append(Interval(each_gene,start,end))
            self.featureIdxs_minus[each_chrom] = IntervalTree(inputlist)

        for each_chrom in temp_nostrand:
            inputlist = []
            for each_gene in temp_nostrand[each_chrom]:
                for (start,end) in temp_nostrand[each_chrom][each_gene]:
                    inputlist.append(Interval(each_gene,start,end))
            self.featureIdxs_nostrand[each_chrom] = IntervalTree(inputlist)


    def getFeatures(self) :
        return self.features

    def Gene_annotation(self,itv_list):
        genes = []
        #fs = []
        for itv in itv_list :
            fs = []
            try:
                if itv[3] == "+" :
                    if itv[0] in self.featureIdxs_plus :
                        fs = self.featureIdxs_plus[itv[0]].find_gene(itv[1],itv[2])


                if itv[3] == "-" :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])


                if itv[3] == "." :
                        if itv[0] in self.featureIdxs_minus:
                            fs = self.featureIdxs_minus[itv[0]].find_gene(itv[1], itv[2])

                        if itv[0] in self.featureIdxs_plus :
                            fs += self.featureIdxs_plus[itv[0]].find_gene(itv[1],itv[2])
                        if itv[0] in self.featureIdxs_nostrand :
                            fs += self.featureIdxs_nostrand[itv[0]].find_gene(itv[1],itv[2])

                if len(fs) > 0:
                        genes = genes + fs

            except:
                    raise


        return genes
