

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

import math
import sys
import array
import operator

from TEToolkit.Constants import *

def normalizeMeans(meansIn):
    total = sum(meansIn)
    meansOut = [0]*len(meansIn)
    sys.stderr.write("total means = " + repr(int(total)) +"\n")
    if total > 0 :
        meansOut = [1.0*x/total for x in meansIn]
   #     for i in range(len(means)) :
   #         means[i] = 1.0*means[i]/total

    #sys.stderr.write("after normalization total means = "+str(sum(meansOut))+"\n")
    return meansOut

def EMUpdate(meansIn, te_features,uniq_counts,multi_reads,estimatedReadLength):
    # reassign multi-reads proportionally to the relative abundance of each TE
    te_counts = {}
    meansOut = [0] *len(meansIn)

    multi_counts = computeAbundances(meansIn,multi_reads)
    sys.stderr.write("total multi counts = "+ repr(int(sum(multi_counts)))+"\n")
    for tid in range(len(meansIn)) :
        tlen = te_features.getLength(tid)
        if tlen <0 :
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
        effectiveLength = tlen - estimatedReadLength + 1
        if effectiveLength > 0 :
            meansOut[tid] = (uniq_counts[tid] + multi_counts[tid])/(1.0*effectiveLength)
        else :
         #   sys.stderr.write("effective length is less than read length\n")
            meansOut[tid] = 0.0


    meansOut = normalizeMeans(meansOut)
    sys.stderr.write("after normalization total means = "+repr(sum(meansOut))+"\n")

    return meansOut

#def sum_(v) :
#    total = 0.0

#    for vid in v :
#        total += vid

#    return total

### ????
def logLikelihood_(means) :

 #   likelihoods = [0.0] *(transcriptsForKmer_.size())

    #Compute the log-likelihood
 #   for kid in transcriptsForKmer :
        #[&likelihoods, &means, this](const BlockedIndexRange& range) ->void {
 #       kmerLikelihood = 0.0
 #       for tid in transcriptsForKmer_[kid] :
 #           kmerLikelihood += transcripts_[tid].binMers[kid] * (means[tid] / transcripts_[tid].length)
  #          if kmerLikelihood > 1e-20 :
  #              likelihoods[kid] = log(kmerLikelihood)

            #else :
                #  " has probability too low\n"
    return 2 #sum_(likelihoods)

def expectedLogLikelihood_(means, counts, te_features) :
    sampProbs = [0.0]*len(means)

    for tid in range(len(counts)) :
        relativeAbundance = means[tid]
   #     sampProbs[tid] = te_features[tid] * relativeAbundance

  #  normalize_(sampProbs)

    if len(means) > 0 :
        return logLikelihood( sampProbs) / len(means)
    else :
        return 0.0

def absdiff_(v0, v1) :
    diff = 0.0
    #for i in range(len(v0)) :
    #    diff += math.fabs(v0[i]-v1[i])
    diff = sum(map(lambda x,y: abs(x-y),v0,v1))
    return diff

def dotProd_(u, v) :
    dot = 0.0
    dot = sum(map(operator.mul,u,v))
    #for i in range(len(u)) :
    #    dot += u[i] * v[i]

    return dot


def EMestimate(te_features,multi_reads,uniq_counts,multi_counts,numItr,estimatedReadLength):

    #estimate average read length
   # estimatedReadLength = averageReadLength(multi_reads)

    #density per base
    sys.stderr.write("multi-reads = %s " % repr(int(len(multi_reads))))
    means0 = []
    for tid in range(len(uniq_counts)) :
        tlen = te_features.getLength(tid)
        if tlen < 0 :
            sys.stderr.write("Error in optimization: the TE does not exist!\n")
            raise
        effectiveLength = tlen  - estimatedReadLength + 1
        #if effectiveLength < 0 :
        #    effectiveLength = 1
        if effectiveLength > 0 :
            means0.append(1.0*(uniq_counts[tid] + multi_counts[tid])/effectiveLength)
        else :
            #sys.stderr.write("effective length is less that read length mean0.\n")
            means0.append(0.0)

    # relative abundance
    means0 = normalizeMeans(means0)
    sys.stderr.write("after normalization total means0 = "+repr(sum(means0))+"\n")
    '''
    /**
         * Defaults for these values taken from the R implementation of
         * [SQUAREM](http://cran.r-project.org/web/packages/SQUAREM/index.html).
         */
    '''
    minStep0 = 1.0
    minStep = 1.0
    maxStep0 = 1.0
    maxStep = 1.0
    mStep = 4.0
    nonMonotonicity = 1.0
    negLogLikelihoodOld = float("inf")
    negLogLikelihoodNew = float("inf")
    # Right now, the # of iterations is fixed, but termination should
    # also be based on tolerance

    #while feval < numItr :
    cur_iter = 0
    #for iter in range(0,numItr) :
    t_size = len(uniq_counts)
    r = [0]*t_size
    r2 = [0] * t_size
    v = [0]*t_size
    meansPrime = [0.0] * t_size
    #means2 = [0.0] * t_size
    #means1 = [0.0] * t_size
    outerIteration = 1
    while cur_iter < numItr :
        cur_iter += 1
        sys.stderr.write("SQUAREM iteraton [" + str(cur_iter) + "]\n")
        sys.stderr.write("1/3\n")
        try :
            means1 = EMUpdate(means0, te_features,uniq_counts,multi_reads,estimatedReadLength)
        except :
            sys.stderr.write("Error in EMupdate\n")
            raise

        if negLogLikelihoodOld != float("inf") :
            negLogLikelihoodOld = -expectedLogLikelihood_(means0)

        sys.stderr.write("2/3\n")
        try:
            means2 = EMUpdate(means1, te_features,uniq_counts,multi_reads,estimatedReadLength)
        except :
            sys.stderr.write("Error in EMupdate\n")
            raise

#        delta = absdiff_(means1, means2)
 #       sys.stderr.write("delta = " + delta + "\n")

        for tid in range(len(means0)) :
            r[tid] = means1[tid] - means0[tid]
            r2[tid] = means2[tid] - means1[tid]
            v[tid] = (means2[tid] - means1[tid]) - r[tid]
        #    sys.stderr.write("mean1 = "+str(means1[tid])+"\t means2= "+str(means2[tid])+"\t"+str(v[tid])+"\n")
        #sys.exit(0)

        rNorm = math.sqrt(dotProd_(r,r))
        r2Norm = math.sqrt(dotProd_(r2,r2))
        vNorm = math.sqrt(dotProd_(v,v))
        rr = dotProd_(r,v)
        rvNorm = math.sqrt(abs(rr))

        if vNorm == 0 :
            means0 = means1
            sys.stderr.write("at iteration " + str(cur_iter) + " vNorm == 0 \n")
            break
        alphaS = rNorm / rvNorm
        alphaS = max(minStep, min(maxStep, alphaS))

        if rNorm < OPT_TOL :
            sys.stderr.write("rNome = OPT_TOL \n")
            break
        if r2Norm < OPT_TOL :
            sys.stderr.write("r2Nome = OPT_TOL \n")
            means0 = means2
            break

        for tid in range(len(means0)) :
            meansPrime[tid] = max(0.0, means0[tid] + 2*alphaS*r[tid] + (alphaS*alphaS)*v[tid])

        # Stabilization step
        if abs(alphaS - 1.0) > 0.01 :
            sys.stderr.write("alpha = " + repr(alphaS) + ".\n ")
            sys.stderr.write("Performing a stabilization step.\n")
            try :
                meansPrime = EMUpdate(meansPrime, te_features,uniq_counts,multi_reads,estimatedReadLength)
            except :
                sys.stderr.write("Error in EMupdate\n")
                raise
             #/** Check for an error in meansPrime **/
             #/** If there is **/
            #if nonMonotonicity == float("inf") :
            #    sys.stderr.write("......in nonMonotonicity 1.......\n")
            #    negLogLikelihoodNew = - 1.0 * expectedLogLikelihood_(meansPrime)
            #else :
            #    sys.stderr.write("......in nonMonotonicity 2.......\n")
            #    negLogLikelihoodNew = negLogLikelihoodOld

            #if negLogLikelihoodNew > negLogLikelihoodOld + nonMonotonicity :
                #exchange values
            #    sys.stderr.write("......in nonMonotonicity 3.......\n")
            #    meansPrime, means2 = means2, meansPrime
            #    negLogLikelihoodNew = - 1.0 * expectedLogLikelihood_(meansPrime)

                if alphaS == maxStep :
                    maxStep = max(maxStep0, maxStep/mStep)
                    alphaS = 1.0

        sys.stderr.write("alpha = " + repr(alphaS) + ", ")

        if alphaS == maxStep :
            maxStep = mStep * maxStep

        if minStep < 0 and alphaS == minStep :
            minStep = mStep * minStep

        meansPrime, means0 = means0, meansPrime

        if not math.isnan(negLogLikelihoodNew) :
            negLogLikelihoodOld = negLogLikelihoodNew

    if cur_iter >= numItr :
        sys.stderr.write("not converge.....\n")
    else :
        sys.stderr.write("converge at iteration " + str(cur_iter)+"\n")
    new_multi_counts = computeAbundances(means0,multi_reads)
    return new_multi_counts

def computeAbundances(meansIn,multi_reads):

    size = len(meansIn)
    multi_counts = [0] * size

    sys.stderr.write("num of multi reads = "+repr(int(len(multi_reads)))+"\n")
    for kid in range(len(multi_reads)) :
        TE_transcripts = multi_reads[kid]
        totalMass = 0.0
        for tid in TE_transcripts :
              totalMass += meansIn[tid]

        if totalMass > 0.0 :
                  norm = 1.0 / totalMass
        else :
                  norm = 0.0

        for tid in TE_transcripts :
              multi_counts[tid] += meansIn[tid] * norm

    sys.stderr.write("total multi counts = "+ repr(int(sum(multi_counts)))+"\n")
    return multi_counts
