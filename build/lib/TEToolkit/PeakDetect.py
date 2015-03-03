

"""Module Description

Copyright (c) 2014, Ying Jin <yjin@cshl.edu >


This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).

@author:  Ying Jin
@contact: yjin@cshl.edu
"""
import os
from math import log as mathlog
from math import ceil, floor
from array import array
import subprocess

from TEToolkit.Prob import poisson_cdf,poisson_cdf_inv
from TEToolkit.Constants import *
from TEToolkit.IO.FeatIO import FWTrackII

class PeakDetect:
    """Class to do the peak calling.

    e.g:
    >>> from MACS14.PeakDetect import PeakDetect
    >>> pd = PeakDetect(treat=treatdata, control=controldata, pvalue=pvalue_cutoff, d=100, scan_window=200, gsize=3000000000)
    >>> pd.call_peaks()
    >>> print pd.toxls()
    """
    def __init__ (self,opt=None,treat=None, control=None,scal_factors=None):
        """Initialize the PeakDetect object.

        """
        self.sf = scal_factors
        self.opt  = opt
        self.info = opt.info
        self.debug = opt.debug
        self.warn = opt.warn

        self.treat = treat
        self.control = control
        self.ratio_treat2control = None
        self.peaks = {}
        self.final_peaks = {}
        self.final_negative_peaks = {}

                
        self.pvalue = opt.log_pvalue
        self.d = opt.fragsize
        self.shift_size = self.d/2
        self.scan_window = opt.scanwindow
        self.gsize = opt.gsize
        
        self.nolambda = False #opt.nolambda

        self.sregion = 1000 #opt.smalllocal
        self.lregion = 10000 #opt.largelocal

        if (self.nolambda):
            self.info("#3 !!!! DYNAMIC LAMBDA IS DISABLED !!!!")

        #self.save_score = opt.store_score
        #self.zwig_tr = opt.zwig_tr
        #self.zwig_ctl= opt.zwig_ctl

    def call_peaks (self):
        """Call peaks function.

        Scan the whole genome for peaks. RESULTS WILL BE SAVED IN
        self.final_peaks and self.final_negative_peaks.
        """
               # w/ control
        self.peaks = self.__call_peaks_w_control ()
       
        return None



    def toxls (self):
        """Save the peak results in a tab-delimited plain text file
        with suffix .xls.
        
        """
        text = "chr\t"
        if self.control and self.peaks:
            text += "\t".join(("start", "end",  "length",  "summit", "tags", "-10*log10(pvalue)", "fold_enrichment", "FDR(%)"))+"\n"
        elif self.peaks:
            text += "\t".join(("start", "end",  "length",  "summit", "tags", "-10*log10(pvalue)", "fold_enrichment"))+"\n"
        else:
            return ""
        
        chrs = self.peaks.keys()
        chrs.sort()
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                text += "%s\t%d\t%d\t%d" % (chrom,peak[0]+1,peak[1],peak[2])
                peak_summit_relative_pos = peak[3]-peak[0]
                text += "\t%d" % (peak_summit_relative_pos)
                text += "\t%d\t%.2f" % (peak[5],peak[6])
                text += "\t%.2f" % (peak[7])
                if self.control:
                    if peak[8]>=100:
                        text += "\t100"
                    else:
                        text += "\t%.2f" % (peak[8])
                text+= "\n"
        return text

 
    def tobed (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n += 1
                text+= "%s\t%d\t%d\tMACS_peak_%d\t%.2f\n" % (chrom,peak[0],peak[1],n,peak[6])
        return text
    
    def towig (self,filename,span=100,libsize=10000000):
        if os.path.isfile(filename) :
            subprocess.call(["rm -f ",filename ],shell=True)
            
        f = open(filename,'w')
        
        f.write("track type=wiggle_0 name= " + self.treat.get_name()+"\n")
        
        for chr in self.peaks.keys() :
            f.write("variableStep chrom="+chr+" span=" + str(span)+"\n") 
            if len(self.peaks[chr]) > 0 : 
                (tags,cnts) = self.treat.get_locations_by_chr_v2(chr)
                peak_idx = 0
                tag_idx = 0
                cur_w = 0
                cur_pos = 0
                while peak_idx < len(self.peaks[chr]) :
                    peak_start = self.peaks[chr][peak_idx][0]
                    peak_end = self.peaks[chr][peak_idx][1]
                    peak_start_bin = peak_start / span
                    peak_end_bin = peak_end/span + 1
                    
                    tag_bin = tags[tag_idx]/span
                    while tag_bin < peak_start_bin and tag_idx < len(tags):
                        tag_idx += 1
                        tag_bin = tags[tag_idx]/span
                    
                    while tag_bin >= peak_start_bin and tag_bin <= peak_end_bin :
                        cur_w += cnts[tag_idx]
                        tag_idx += 1
                        tag_bin = tags[tag_idx]/span
                    
                    if cur_w > 0 :
                        cur_w = cur_w * self.sf[0]*libsize/self.treat.libsize()
                        f.write(str(peak_start_bin)+"\t"+str(cur_w)+"\n")
                        cur_w = 0
                        
                    peak_idx += 1                
                
        f.close()
        
        return 

    def get_counts_by_chr(self,chrom,strand) :
        
        return self.treat.get_counts_by_chr(chrom,0)
    
    def clean(self, chr):
        return self.treat.clean(chr)
    
    def summitsToBED (self):
        text = ""
        chrs = self.peaks.keys()
        chrs.sort()
        n = 0
        for chrom in chrs:
            for peak in self.peaks[chrom]:
                n += 1
                text+= "%s\t%d\t%d\tMACS_peak_%d\t%.2f\n" % (chrom,peak[3]-1,peak[3],n,peak[4])
        return text

    def __add_fdr (self, final, negative): 
        """
        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type
        """
        pvalue2fdr = {}
        pvalues_final = []
        pvalues_negative = []
        chrs = final.keys()
        a = pvalues_final.append
        for chrom in chrs:
            for i in final[chrom]:
                a(i[6]) # i[6] is pvalue in peak info
                pvalue2fdr[i[6]]=None
        chrs = negative.keys()
        a = pvalues_negative.append
        for chrom in chrs:
            for i in negative[chrom]:
                a(i[6])
        pvalues_final.sort(reverse=True)
        pvalues_final_l = len(pvalues_final)
        pvalues_negative.sort(reverse=True)
        pvalues_negative_l = len(pvalues_negative)        
        pvalues = pvalue2fdr.keys()
        pvalues.sort(reverse=True)
        index_p2f_pos = 0
        index_p2f_neg = 0
        for p in pvalues:
            while index_p2f_pos<pvalues_final_l and p<=pvalues_final[index_p2f_pos]:
                index_p2f_pos += 1
            n_final = index_p2f_pos

            while  index_p2f_neg<pvalues_negative_l and p<=pvalues_negative[index_p2f_neg]:
                index_p2f_neg += 1
            n_negative = index_p2f_neg
            pvalue2fdr[p] = 100.0 * n_negative / n_final

        new_info = {}
        chrs = final.keys()
        for chrom in chrs:
            new_info[chrom] = []
            for i in final[chrom]:
                tmp = list(i)
                tmp.append(pvalue2fdr[i[6]])
                new_info[chrom].append(tuple(tmp))      # i[6] is pvalue in peak info
        return new_info

    def __call_peaks_w_control (self):
        """To call peaks with control data.

        A peak info type is a: dictionary

        key value: chromosome

        items: (peak start,peak end, peak length, peak summit, peak
        height, number of tags in peak region, peak pvalue, peak
        fold_enrichment) <-- tuple type
        """
        #self.lambda_bg = float(self.scan_window)*self.treat.total/self.gsize
        #self.debug("#3 background lambda: %.2f " % (self.lambda_bg))
        #self.min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),self.lambda_bg)+1
        #self.debug("#3 min tags: %d" % (self.min_tags))

        self.ratio_treat2control = float(self.treat.total)/self.control.total
        if self.ratio_treat2control > 2 or self.ratio_treat2control < 0.5:
            self.warn("Treatment tags and Control tags are uneven! FDR may be wrong!")
        self.info("#3 shift treatment data")
        self.__shift_trackI(self.treat)
        self.info("#3 merge +/- strand of treatment data")

        self.treat.merge_plus_minus_locations_naive ()

        self.info("#3 after shift and merging, tags: %d" % (self.treat.total))
 #       self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #if self.opt.wigextend:
            #    zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.single_profile)
            #else:
   #     zwig_write(self.treat,self.opt.wig_dir_tr,self.zwig_tr,self.d,log=self.info,space=self.opt.space,single=self.opt.single_profile)
        self.info("#3 call peak candidates")
        #peak_candidates = self.__call_peaks_from_trackI (self.treat)

        peak_candidates = self.__call_peaks_from_trackI_v2 (self.treat) #Ying
        
        self.info("#3 shift control data")
        self.info("#3 merge +/- strand of control data")
        self.__shift_trackI(self.control)
        self.control.merge_plus_minus_locations_naive ()

        self.info("#3 after shift and merging, tags: %d" % (self.control.total))
  #      self.info("#3 save the shifted and merged tag counts into wiggle file...")
            #if self.opt.wigextend:
            #    zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.opt.wigextend,log=self.info,space=self.opt.space,single=self.opt.single_profile)
            #else:
   #     zwig_write(self.control,self.opt.wig_dir_ctl,self.zwig_ctl,self.d,log=self.info,space=self.opt.space,single=self.opt.single_profile)
        self.info("#3 call negative peak candidates")
#        negative_peak_candidates = self.__call_peaks_from_trackI (self.control)

        negative_peak_candidates = self.__call_peaks_from_trackI_v2 (self.control) # Ying

        # build score
        #if self.save_score:
        #    self.info("#4 build scores")
        #    scoreswig = self.__build_score_wigtrackI(treatwig,controlbkI,self.d,space=self.opt.space,bglambda=self.lambda_bg)
        #    zwigfile = file(self.opt.name+".score.wig","w")
        #    self.info("#4.1 save scores to wiggle file")            
        #    scoreswig.write_wig(zwigfile,"score")
        #    self.info("compress the wiggle file using gzip...")
        #    os.system("gzip "+self.opt.name+".score.wig")
        
        self.info("#3 use control data to filter peak candidates...")
        #self.final_peaks = self.__filter_w_control(peak_candidates,self.treat,self.control, self.ratio_treat2control,fake_when_missing=True)
        sf = [1,1]
        sf[0] = self.sf[0]
        sf[1] = self.sf[1]
        self.final_peaks = self.__filter_w_control_v2(sf,peak_candidates,self.treat,self.control, fake_when_missing=True, to_small_sample=False)
        self.info("#3 find negative peaks by swapping treat and control")

        #self.final_negative_peaks = self.__filter_w_control(negative_peak_candidates,self.control,self.treat, 1.0/self.ratio_treat2control,fake_when_missing=True)
        sf = [1,1]
        sf[0] = self.sf[1]
        sf[1] = self.sf[0]
        self.final_negative_peaks = self.__filter_w_control_v2(sf,negative_peak_candidates,self.control,self.treat, fake_when_missing=True, to_small_sample=False)
        return self.__add_fdr (self.final_peaks, self.final_negative_peaks)

  


    def __filter_w_control_v2 (self, sf,peak_info, treatment, control, pass_sregion=False, write2wig= False, fake_when_missing=False, to_small_sample=False ):
        """Use control data to calculate several lambda values around
        1k, 5k and 10k region around peak summit. Choose the highest
        one as local lambda, then calculate p-value in poisson
        distribution.

        Parameters:

        1. pass_sregion: If set True, the slocal lambda will be
        ignored. Use this when the control is not available.
        
        2. write2wig: obselete
        
        3. fake_when_missing: when a chromosome is missing in control
        but existing in IP or vice versa, MACS will fake a tag to pass
        the process.
        
        4. to_small_sample: when set as True, balance the number of
        tags by linearly scaling larger sample to smaller sample. The
        default behaviour is to linearly scale smaller to larger one.

        Return value type in this format:
        a dictionary
        key value : chromosome
        items : array of (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags,peak_pvalue,peak_fold_enrichment)
        """
        lambda_bg0 = float(self.scan_window)*treatment.total/self.gsize # bug fixed...
        
 #       if treatment.total>control.total:
 #           t_ratio = 1.00
 #           c_ratio = float(treatment.total)/control.total
 #       else:
 #           t_ratio = float(control.total)/treatment.total
 #           c_ratio = 1.00

 #       if to_small_sample:
 #           tmp = t_ratio
 #           t_ratio = 1/c_ratio
 #           c_ratio = 1/tmp
        
        t_ratio = sf[0]
        c_ratio = sf[1]
        
   #     self.info("t_ratio %s" % (t_ratio))
   #     self.info("c_ratio %s" % (c_ratio))
        
        final_peak_info = {}
        chrs = peak_info.keys()
        chrs.sort()
        total = 0
        for chrom in chrs:
          #  self.info("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            final_peak_info[chrom] = []
            peak_list = peak_info[chrom]
            try:
                (ctags,ccnts) = control.get_locations_by_chr_v2(chrom)
                #ccnts = control.get_counts_by_chr(chrom,0)  
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ctags = [-1,]
                    ccnts = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
            try:
                (ttags,tcnts) = treatment.get_locations_by_chr_v2(chrom)
                #tcnts = treatment.get_counts_by_chr(chrom,0)
            except:
                self.warn("Missing %s data, skip it..." % (chrom))
                if fake_when_missing:
                    ttags = [-1,]
                    tcnts = [-1,]
                    self.warn("Fake a tag at %s:%d" % (chrom,-1))
                    tmp=[]
                else:
                    continue
         #   self.info("ttags size %d" % (len(ttags)))  
         #   self.info("ctags size %d" % (len(ttags)))  
         #   self.info("ttags size %d" % (len(tcnts)))  
         #   self.info("ctags size %d" % (len(ccnts)))  
            
            index_ctag = 0      # index for control tags
            index_ttag = 0      # index for treatment tags
            flag_find_ctag_locally = False
            flag_find_ttag_locally = False            
            prev_index_ctag = 0
            prev_index_ttag = 0            
            len_ctags =len(ctags)
            len_ttags =len(ttags)            
            for i in range(len(peak_list)):
                (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags) = peak_list[i]

                #window_size_4_lambda = min(self.first_lambda_region,max(peak_length,self.scan_window))
                window_size_4_lambda = max(peak_length,self.scan_window)
                lambda_bg = lambda_bg0/self.scan_window*window_size_4_lambda*t_ratio               
                if self.nolambda:
                    # skip local lambda
                    local_lambda = lambda_bg
                    tlambda_peak = float(peak_num_tags)/peak_length*window_size_4_lambda
                else:
                    left_peak = peak_start+self.shift_size # go to middle point of the first fragment
                    right_peak = peak_end-self.shift_size  # go to middle point of the last fragment
                    left_lregion = peak_summit-self.lregion/2
                    left_sregion = peak_summit-self.sregion/2
                    right_lregion = peak_summit+self.lregion/2
                    right_sregion = peak_summit+self.sregion/2
                    #(cnum_10k,cnum_5k,cnum_1k,cnum_peak) = (0,0,0,0)
                    #(tnum_10k,tnum_5k,tnum_1k,tnum_peak) = (0,0,0,0)
                    (cnum_sregion, cnum_lregion, cnum_peak, tnum_sregion, tnum_lregion, tnum_peak) = (0,0,0,0,0,0)
                    #smallest = min(left_peak,left_10k,left_5k,left_1k)
                    #largest = max(right_peak,right_10k,right_5k,right_1k)

                    while index_ctag < len_ctags:
                        if ctags[index_ctag] < left_lregion:
                            # go to next control tag
                            index_ctag+=1
                        elif index_ctag+1 >= len_ctags or right_lregion < ctags[index_ctag]:
                            # If move outof the lregion or reach the chromosome end
                            # finalize and go to next peak region
                            # Thanks to Jake Biesinger
                            flag_find_ctag_locally = False
                            index_ctag = prev_index_ctag 
                            break
                        else:
                            if not flag_find_ctag_locally:
                                flag_find_ctag_locally = True
                                prev_index_ctag = index_ctag
                            p = ctags[index_ctag]
                            c = ccnts[index_ctag]
                            if left_peak <= p <= right_peak:
                                cnum_peak += c
                            if left_sregion <= p <= right_sregion:
                                cnum_sregion +=c
                                cnum_lregion +=c
                            else:
                                cnum_lregion += c
                            index_ctag += 1 # go to next tag

                    while index_ttag < len_ttags:
                        if ttags[index_ttag] < left_lregion:
                            # go to next treatment tag
                            index_ttag+=1
                        elif index_ttag+1 >= len_ttags or right_lregion < ttags[index_ttag]:
                            # If move outof the lregion or reach the chromosome end
                            # finalize and go to next peak region
                            # Thanks to Jake Biesinger
                            flag_find_ttag_locally = False
                            index_ttag = prev_index_ttag 
                            break
                        else:
                            if not flag_find_ttag_locally:
                                flag_find_ttag_locally = True
                                prev_index_ttag = index_ttag
                            p = ttags[index_ttag]
                            c = tcnts[index_ttag]
                            if left_peak <= p <= right_peak:
                                tnum_peak +=c
                            if left_sregion <= p <= right_sregion:
                                tnum_sregion +=c
                                tnum_lregion += c
                            else:
                                tnum_lregion += c
                            index_ttag += 1 # go to next tag

                    clambda_peak = float(cnum_peak)/peak_length*c_ratio*window_size_4_lambda

                    clambda_lregion = float(cnum_lregion)/self.lregion*c_ratio*window_size_4_lambda

                    clambda_sregion = float(cnum_sregion)/self.sregion*c_ratio*window_size_4_lambda

                    tlambda_peak = float(tnum_peak)/peak_length*t_ratio*window_size_4_lambda

                    tlambda_lregion = float(tnum_lregion)/self.lregion*t_ratio*window_size_4_lambda

                    tlambda_sregion = float(tnum_sregion)/self.sregion*t_ratio*window_size_4_lambda

                    if pass_sregion:
                        # for experiment w/o control, peak region lambda and sregion region lambda are ignored!
                        local_lambda = max(lambda_bg,tlambda_lregion)
                    else:
                        # for experiment w/ control
                        local_lambda = max(lambda_bg,clambda_peak,clambda_lregion,clambda_sregion)

                #print(local_lambda)
               # if local_lambda == 0 :
               #     local_lambda = 0.001
                p_tmp = poisson_cdf(tlambda_peak,local_lambda,lower=False)
                if p_tmp <= 0:
                    peak_pvalue = 3100
                else:
                    peak_pvalue = mathlog(p_tmp,10) * -10

                if peak_pvalue > self.pvalue:
                    n_chrom += 1
                    total += 1
                    peak_fold_enrichment = float(peak_height)/local_lambda*window_size_4_lambda/self.d
                    final_peak_info[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags,peak_pvalue,peak_fold_enrichment))
                # uncomment the following two lines, MACS will report the peaks been rejected.    
                #else:
                #    #self.debug("Reject the peak at %s:%d-%d with local_lambda: %.2f and -log10pvalue: %.2f" % (chrom,peak_start,peak_end,local_lambda,peak_pvalue))

            self.debug("#3 peaks whose pvalue < cutoff: %d" % (n_chrom))
        self.info("#3 Finally, %d peaks are called!" % (total))
        return final_peak_info
    
    
    def __call_peaks_from_trackI_v2 (self, trackI):
        """ Ying Call peak candidates from trackI data. Using every tag as
        step and scan the self.scan_window region around the tag. If
        tag number is greater than self.min_tags, then the position is
        recorded.

        Return: data in this format. (peak_start,peak_end,peak_length,peak_summit,peak_height,peak_num_tags)
        """
        lambda_bg0 = float(self.scan_window)*trackI.total/self.gsize # bug fixed...        
        self.info("#3 all peaks candidate lambda_bg0 = : %f" % (lambda_bg0))       
        
        min_tags = poisson_cdf_inv(1-pow(10,self.pvalue/-10),lambda_bg0)+1
        peak_candidates = {}
    
        self.info("#3 search peak condidates... min_tags %d" % (min_tags))
        #self.info("#3 search peak condidates... scan window %d" % (self.scan_window))
        chrs = trackI.get_chr_names()
        total = 0
        for chrom in chrs:
        #    self.info("#3 Chromosome %s" % (chrom))
            n_chrom = 0
            peak_candidates[chrom] = []
            (tags,cnts) = trackI.get_locations_by_chr_v2(chrom)
          #  for k in range(len(tags)) :
          #      print(str(tags[k])+"\t"+str(cnts[k]))
            #cnts = trackI.get_counts_by_chr(chrom,0)

            len_t = len(tags)
            if len_t ==0 :
                continue
 
            cpr_tags = []       # Candidate Peak Region tags
            #cpr_tags.extend(tags[:min_tags-1])
            cpr_cnts = []     # Candidate peak region counts
            i = -1
            number_cpr_tags = 0
            while number_cpr_tags <= min_tags-1:
                i += 1
                if(i >= len(cnts)):
                    break
                number_cpr_tags += cnts[i]
                cpr_cnts.append(cnts[i])
                cpr_tags.append(tags[i])
                #i += 1
            #p = min_tags-1 # Next Tag Index
            if i == len(cnts) :
                continue
            if len(cpr_tags) > 0 :
                cpr_tags.pop(i)
                cpr_cnts.pop(i)
                number_cpr_tags -= cnts[i]
            p = i
            cur_num = i 

            while p < len_t:

                if number_cpr_tags >= min_tags:
                 
                    if tags[p] - cpr_tags[-1*cur_num+1] <= self.scan_window:
                        # add next tag, if the new tag is less than self.scan_window away from previous no. min_tags tag
                        cpr_tags.append(tags[p])
                        cpr_cnts.append(cnts[p])
                        number_cpr_tags += cnts[p]
                        #cur_num +=1
                        p+=1
                    else:
                        # candidate peak region is ready, call peak...
			#self.info("detect a candidate peak")
        
                        (peak_start,peak_end,peak_length,peak_summit,peak_height) = self.__tags_call_peak_v2 (cpr_tags,cpr_cnts)
                        peak_candidates[chrom].append((peak_start,peak_end,peak_length,peak_summit,peak_height,round(number_cpr_tags)))

                        cpr_tags = [tags[p]] # reset
                        cpr_cnts = [cnts[p]]
                        number_cpr_tags = cnts[p]
                        cur_num = 1
                        #last_idx = p
                        #last_cnt = cnts[p]
                        p += 1               
                        
                       
                        total += 1
                        n_chrom += 1

                else:
                    # add next tag, but if the first one in cpr_tags
                    # is more than self.scan_window away from this new
                    # tag, remove the first one from the list
                    if tags[p] - cpr_tags[0] >= self.scan_window:
                        number_cpr_tags -= cpr_cnts[0]
                        cpr_tags.pop(0)
                        cpr_cnts.pop(0)
                        cur_num -= 1
                        
                    cpr_tags.append(tags[p])
                    cpr_cnts.append(cnts[p])
                    number_cpr_tags += cnts[p]
                    cur_num += 1
                    #last_idx +=1
                    #last_cnt = cnts[last_idx]
                    p+=1
          #  self.info("#3 peak candidates: %d" % (n_chrom))
        self.info("#3 Total number of candidates: %d" % (total))
        return self.__remove_overlapping_peaks(peak_candidates)
                
            
    def __tags_call_peak_v2 (self, tags,tags_cnts ):
        """Project tags to a line. Then find the highest point.

        """
        start =  tags[0]-self.d/2
        end = tags[-1]+self.d/2       # +1 or not?
        region_length = int(end - start)
        line= [0]*region_length
        for k in range(len(tags)):
            tag = int(tags[k])
            tag_projected_start = tag-start-self.d/2
            tag_projected_end = tag-start+self.d/2
            for i in range(tag_projected_start,tag_projected_end):
                line[i] += tags_cnts[k]
        tops = []
        top_height = 0
        for i in range(len(line)):
            if line[i] > top_height:
                top_height = line[i]
                tops = [i]
            elif line[i] == top_height:
                tops.append(i)
        peak_summit = int(tops[len(tops)/2]+start)
        return (start,end,region_length,peak_summit,top_height)
    
    def __shift_trackI (self, trackI):
        """Shift trackI data to right (for plus strand) or left (for
        minus strand).

        """
        chrs = trackI.get_chr_names()
        number_removed_tags = 0
        for chrom in chrs:
            tags = trackI.get_locations_by_chr(chrom)
           # cnts = trackI.get_counts_by_chr(chrom,0)
           # self.info("shift chromosome %s" % (chrom))
            #self.info("length of chromosome strand plus: %d" % (len(tags[0])))
            #self.info("counts of chromosome strand plus: %d" % (len(cnts)))
	    # plus
            for i in range(len(tags[0])):
                tags[0][i]+=self.shift_size
                #cnts[i].pos += self.shift_size
            # minus
	    #cnts = []
            #cnts = trackI.get_counts_by_chr(chrom,1)
            #self.info("length of chromosome strand minus: %d " % (len(tags[1])))
            #self.info("length of chromosome strand minus: %d " % (len(cnts)))
            for i in range(len(tags[1])):
                tags[1][i]-=self.shift_size
                #cnts[i].pos -= self.shift_size
            # remove the tags extended outside of chromosome start
           # self.info("remove tags outside the range ")
            while tags[1]:
                #self.info("tags length : %d " % (len(tags[1])))
                if tags[1][0]-self.shift_size<0:
                    number_removed_tags += 1
		    #cnts.pop(0)
                    tags[1].pop(0)
                else:
                    break

        self.debug("# %d tag(s) extended outside of chromosome start are removed!" % number_removed_tags)
        return

    def __remove_overlapping_peaks (self, peaks ):
        """peak_candidates[chrom] = [(peak_start,peak_end,peak_length,peak_summit,peak_height,number_cpr_tags)...]

        """
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            n_append = new_peaks[chrom].append
            prev_peak = None
            peaks_chr = peaks[chrom]
            for i in xrange(len(peaks_chr)):
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
        return new_peaks
    
 
