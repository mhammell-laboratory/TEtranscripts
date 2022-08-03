'''
Created on Oct 17, 2011

@author: yjin
'''

def base_m(m_ref1,m_ref2):
    target = []

    for i in range(len(m_ref1)):
        target.append((m_ref1[i]+m_ref2[i])/2)

    return target

def raw_count_mean(m_ref1,m_ref2,norm_ref) :
    target_ref = []
    for i in range(len(m_ref1)):
        target_ref.append(raw_sum_mean(m_ref1[i],m_ref2[i],norm_ref))

    return target_ref

def raw_sum_mean(m1,m2,tr_norm):
    v = 0;
    for i in tr_norm :
        v += 1.0/i;

    return v*(m1+m2)/2



def est_var(tr_fit,base_mean,norm,tr_raw_sum_mean):
    #including raw variance scv adjust
    s = 0
    target = []
    for n in norm :
        s += (1/n)**2

    for k in range(len(tr_fit)):
        var = tr_fit[k] * base_mean[k] * base_mean[k] * s + tr_raw_sum_mean[k]
        correct = 1.00000001 * tr_raw_sum_mean[k]

        if var < correct :
            target.append(correct)
        else:
            target.append(var)

    return target
