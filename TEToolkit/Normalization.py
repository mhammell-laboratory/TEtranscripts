'''
Created on Oct 13, 2011

two normalization methods are supported.
sequence depth
bin correlation

@author: Ying Jin
'''
#from rpy2.robjects.packages import importr
#from rpy2.robjects import FloatVector
#import rpy2.robjects as robjects
#import rpy2.robjects.lib.ggplot2 as ggplot2
import subprocess


from array import array
#import numpy as np

import logging
import sys
from TEToolkit.Constants import *
from TEToolkit.ShortRead.ParseBEDFile import *


#stats = importr('stats')
#grdevices = importr('grDevices')


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

def __binCorr2r(filename):
    rfhd=open(filename,'w')

    rfhd.write("#!/bin/env Rscript\n")

    rfhd.write("args = commandArgs(TRUE)\n")

    rfhd.write("infile = args[1]\n")
    rfhd.write("outfile = args[2]\n")

    rfhd.write("d<-read.delim(infile,header=T,stringsAsFactors=F) \n")

    rfhd.write("colnum <- length(d) \n")

    rfhd.write("min_sf_idx = 1 \n")
    rfhd.write("min_sf = 1.0 \n")

    rfhd.write("for (i in 2:colnum) {\n")

    rfhd.write("lm.r<-lm(d[,1]~d[,i]-1)\n")

    rfhd.write("cur_sf <- as.numeric(lm.r$coefficients[1])\n")

    #rfhd.write("if (min_sf > cur_sf ) { \n")
    rfhd.write("min_sf = cur_sf;\n")
    #rfhd.write("min_sf_idx = i }\n")

    rfhd.write("png(outfile,height=5,width=5,res=500,units='in')\n")
    rfhd.write("plot(d[,1]~d[,2],xlab=colnames(d)[1],ylab=colnames(d)[2],log='xy') \n")
    rfhd.write("dev.off() }\n")
    rfhd.write("cat(min_sf)\n")
    rfhd.close()




def __output_count_tbl(list1,list2, smp1,smp2, fname):


    try:
        f = open(fname, 'w')
    except IOError:
        error("Cannot create report file %s !\n" % (fname))
        sys.exit(1)
    else:
        cnt_tbl = dict()
        header = ""+smp1+"\t"+smp2
        f.write(header+"\n")
        for i in range(len(list1)) :
            f.write(str(list1[i]) + "\t" +str(list2[i])+"\n")

        f.close()

    return

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

    #call linear regression function
 #   refsmp = FloatVector(reads[0])
 #   robjects.globalenv['A'] = refsmp

    min_sf_idx = 0
    min_sf = 1.0

    __binCorr2r("bin_corr.r")

    for i in range(tsize+tisize):
        tmpfname = "."+fnames[0]+"-"+fnames[i]
        __output_count_tbl(reads[0],reads[i],fnames[0],fnames[i],tmpfname)
        outfname = prj_name+"_"+fnames[0]+"vs"+fnames[i]+".png"
        if i != 0 :
            #msf = subprocess.check_output(["Rscript","bin_corr.r",tmpfname,outfname])
            msf = subprocess.Popen(["Rscript","bin_corr.r",tmpfname,outfname],stdout=subprocess.PIPE)
            msf = float(msf.communicate()[0])
        else :
            msf = 1
     #   subprocess.call(["rm","-f",tmpfname])
        if min_sf > msf :
            min_sf = msf
            min_sf_idx = i
        sf1.append(msf)

    for i in range(tsize+tisize,csize + cisize + tsize+tisize):
        tmpfname = "."+fnames[0]+"-"+fnames[i]
        outfname = prj_name+"_"+fnames[0]+"vs"+fnames[i]+".png"
        __output_count_tbl(reads[0],reads[i],fnames[0],fnames[i],tmpfname)

        #msf = subprocess.check_output(["Rscript","bin_corr.r",tmpfname,outfname])
        msf = subprocess.Popen(["Rscript","bin_corr.r",tmpfname,outfname],stdout=subprocess.PIPE)
        msf = float(msf.communicate()[0])
      #  subprocess.call(["rm","-f",tmpfname])
        if min_sf > msf :
            min_sf = msf
            min_sf_idx = i
        sf2.append( msf)


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

'''
def plot_sf(t_reads,reads,sel_idx,smps,ref_idx):

    sf1 = []

    # plotting code here
    refsmp = FloatVector(reads[ref_idx])
    robjects.globalenv['A'] = refsmp

    for i in range(len(reads)):
        robjects.globalenv['B'] = FloatVector(reads[i])
        res = stats.lm('A ~ B - 1')
        slm = robjects.r.summary(res)
        sf1.append(res[0][0])
        r_sqr = slm[8][0]

        if i != ref_idx :
            fname = smps[ref_idx]+"_vs_"+ smps[i] + ".png"
            grdevices.png(file=fname, width=512, height=512)

            xlab_desc = 'Reads per 10kbp bin of sample ' + smps[i]
            ylab_desc = 'Reads per 10kbp bin of sample ' + smps[ref_idx]
            tt = "Slope: " + str(round(res[0][0],2)) + "\tR-squared: " + str(round(r_sqr,2))

            d = {'ref': robjects.FloatVector(t_reads[ref_idx]),
            'smp': robjects.FloatVector(t_reads[i]),'background':robjects.IntVector(sel_idx)}
            dataf = robjects.DataFrame(d)

            gp = ggplot2.ggplot(dataf)
            pp = gp + \
                 ggplot2.aes_string(y = 'ref', x = 'smp', col='factor(background)') + \
                 ggplot2.geom_point(size=2,shape=16) + \
                 ggplot2.stat_abline(intercept=0,slope=res[0][0],col='blue') + \
                 ggplot2.opts(title=tt) + \
                 ggplot2.scale_colour_discrete(name="",breaks = [0,1], labels = ['all','selected']) + \
                 ggplot2.scale_x_continuous(xlab_desc) + \
                 ggplot2.scale_y_continuous(ylab_desc)

  #          pp.plot()
 #           grdevices.dev_off()

#    return sf1
'''
