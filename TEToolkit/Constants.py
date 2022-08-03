'''
Created on Oct 12, 2011

@author: Ying Jin
'''
WIN_SIZE = 1000
GAP = 1000
PRJ_NAME = 'TEpeaks_out'
STEP = 100        # default step/bin size.
FRAG_SIZE = 200    # default fragment size.
#$range = 100;
CUTOFF = 1
epsilon = 0.00000001
NN = 0.99
SHIFTSIZE = 100
MAX_PAIRNUM = 1000

HS_CHROM = 'data/hg19_chromsize'
RN_CHROM = 'data/rn4_chromsize'
MM_CHROM = 'data/mm9_chromsize'
DM_CHROM = 'data/dm3_chromsize'
MAX_BIT = 50000

TEindex_BINSIZE = 500
OPT_TOL = 0.0001 #tolerance for iterative optimization

efgsize = {"hg":2.7e9,
           "mm":1.87e9,
           "tm":0.789e9,
           "dm":1.2e8}


NORM_METHOD = 'sd'
STAT_METHOD = 'gt'
BIN_SIZE = 10000 # for computing bin correlation
P_VAL = 1e-5


hg19_chrom_lengths = {'chr1':249250621,  'chr2':243199373, 'chr3':198022430,
                      'chr4':191154276, 'chr5':180915260, 'chr6':171115067,
                      'chr7':159138663,  'chr8':146364022, 'chr9':141213431,
                      'chr10':135534747, 'chr11':135006516, 'chr12':133851895,
                      'chr13':115169878, 'chr14':107349540, 'chr15':102531392,
                      'chr16':90354753,  'chr17':81195210,  'chr18':78077248,
                      'chr19':59128983,  'chr20':63025520,  'chr21':48129895,
                      'chr22':51304566,  'chrX':155270560,  'chrY':59373566,
                      'chrM':16571}

mm9_chrom_lengths = {'chr1_random':1231697,'chr3_random':41899,'chr13_random':400311,'chr17_random':628739,'chr4_random':160594,'chrX_random':1785075,
                     'chr5_random':357350,'chr7_random':362490,'chr8_random':849593,'chr9_random':449403,'chrUn_random':5900358,'chrY_random':58682461,
                     'chr1':197195432, 'chr2':181748087, 'chr3':159599783,
                     'chr4':155630120, 'chr5':152537259, 'chr6':149517037,
                     'chr7':152524553, 'chr8':131738871, 'chr9':124076172,
                     'chr10':129993255, 'chr11':121843856, 'chr12':121257530,
                     'chr13':120284312, 'chr14':125194864, 'chr15':103494974,
                     'chr16':98319150, 'chr17':95272651, 'chr18':90772031,
                     'chr19':61342430, 'chrX':166650296, 'chrY':15902555,
                     'chrM':16299}

dm3_chrom_lengths = {'chr2L':23011544,
                        'chr2LHet':368872,
                        'chr2R':21146708,
                        'chr2RHet':3288761,
                        'chr3L':24543557,
                        'chr3LHet':2555491,
                        'chr3R':27905053,
                        'chr3RHet':2517507,
                        'chr4':1351857,
                        'chrX':22422827,
                        'chrXHet':204112,
                        'chrYHet':347038,
                        'chrU':10049037,
                        'chrUextra':29004656,
                        'chrM':19517,
                        'X-TAS':9872}
tm24_chrom_lengths = {'ch00':21805821,
                         'ch01':90304244,
                         'ch02':49918294,
                         'ch03':64840714,
                         'ch04':64064312,
                         'ch05':65021438,
                         'ch06':46041636,
                         'ch07':65268621,
                         'ch08':63032657,
                         'ch09':67662091,
                         'ch10':64834305,
                         'ch11':53386025,
                         'ch12':65486253}

species_chrom_lengths={
                       'mm9':mm9_chrom_lengths,
                       'hg19':hg19_chrom_lengths,
                       'dm3':dm3_chrom_lengths,
                       'tm24':tm24_chrom_lengths};


