#!/usr/bin/env python
"""reportLamp.py converts LAMP result to report.

@author:     ohta
@copyright:  2015 Amelieff Corporation. All rights reserved.
@license:
@contact:
@deffield    updated: Updated

"""
import re
import sys
from optparse import OptionParser


RANK_THRESHOLD_DEFAULT = -sys.maxint
DEFAULT_VALUE = -sys.maxint
SEP = '\t'


def readFiles(lampfile, distfile, expfile, rank_threshold):
    """Read peak files and set flag & distance to gene2peaks.

    Keyword arguments:
    lampfile -- Output of LAMP
    distfile -- Output of checkPeak.py
    expfile  -- Output of checkExp.py
    rank_threshold -- Rank threshold
    Returns: Dictionary, Dictionary, List, List

    """

    gene2dist = {}
    combarr = []
    gene2comb = {}
    r_file = open(lampfile, 'r')
    count = 0
    ecount = 0
    retsu_part = re.compile('^[-+]?\d+(\.\d+)?([eE][-+]?[0-9]+)?$')
    for fh in r_file:
        count += 1
        if fh in ('\n', '\r'):
            ecount += 1
            continue
        if (0 == fh.find('#')) or (0 == fh.find('Rank')) or \
           (0 == fh.find('Time')):
            ecount += 1
            continue

        arr = fh.split('\t')
        if 4 > len(arr):
            print >> sys.stderr, 'Error: Less columns at line ' + \
                str(count) + ' in ' + lampfile
            sys.exit()

        if None is retsu_part.search(arr[0]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 1 in ' + lampfile
            sys.exit()

        if None is retsu_part.search(arr[1]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 2 in ' + lampfile
            sys.exit()

        rank = int(arr[0])
        comb = arr[3]
        labelarr = []

        if rank_threshold == RANK_THRESHOLD_DEFAULT or rank <= rank_threshold:
            combarr.append(comb)

    r_file.close()

    if count - ecount <= 0:
        print >> sys.stderr, 'Error: No valid line in ' + lampfile
        sys.exit()

    count = 0
    ecount = 0
    retsuval = 2
    r2_file = open(distfile, 'r')

    for fh_2 in r2_file:
        count += 1
        if fh_2 in ('\n', '\r'):
            ecount += 1
            continue

        arr_2 = fh_2.split(',')
        if retsuval > len(arr_2):
            print >> sys.stderr, 'Error: Less columns at line ' + \
                str(count) + ' in ' + distfile
            sys.exit()

        gene = arr_2[0]
        arr_2 = arr_2[1:]

        if gene == '#gene':
            retsuval = len(arr_2) + 1
            ecount += 1
            labelarr = arr_2
        else:
            gene2dist[gene] = SEP.join(map(str, arr_2))

            if len(labelarr) != len(arr_2):
                print >> sys.stderr, 'Error: Less columns at line ' + \
                    str(count) + ' in ' + distfile
                sys.exit()

            for comb_2 in combarr:
                if gene not in gene2comb:
                    gene2comb[gene] = {}
                if comb_2 not in gene2comb[gene]:
                    gene2comb[gene][comb_2] = DEFAULT_VALUE

                tfCnt = 0
                tfArr = comb_2.split(',')
                for tf in tfArr:
                    for i in range(len(arr_2)):
                        label = labelarr[i]
                        value = arr_2[i]
                        if tf == label and value != '-':
                          tfCnt = tfCnt + 1

                if tfCnt == len(tfArr):
                  gene2comb[gene][comb_2] = 0

    r2_file.close()

    if count - ecount <= 0:
        print >> sys.stderr, 'Error: No valid line in ' + distfile
        sys.exit()
    count = 0
    ecount = 0
    r3_file = open(expfile, 'r')

    for fh_3 in r3_file:
        count += 1

        if fh_3 in ('\n', '\r'):
            ecount += 1
            continue

        if 0 == fh_3.find('#gene'):
            ecount += 1
            continue

        arr_3 = fh_3.split(',')

        if 2 > len(arr_3):
            print >> sys.stderr, 'Error: Less columns at line ' + \
                str(count) + ' in ' + expfile
            sys.exit()

        if None is retsu_part.search(arr_3[1]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 2 in ' + expfile
            sys.exit()

        gene = str(arr_3[0])
        exp  = int(arr_3[1])
        for comb_3 in combarr:
            if gene2comb[gene][comb_3] != DEFAULT_VALUE and exp > 0:
                gene2comb[gene][comb_3] = 1
    r3_file.close()

    if count - ecount <= 0:
        print >> sys.stderr, 'Error: No valid line in ' + expfile
        sys.exit()

    for i_2 in range(len(labelarr)):
        st = ''.join(map(str, labelarr[i_2]))
        st = st.rstrip()
        labelarr[i_2] = 'Distance_' + ''.join(map(str, st))
    return gene2comb, gene2dist, labelarr, combarr


def main():
    lampFile = ''
    distFile = ''
    expFile = ''
    outFile = ''
    try:
        parser = OptionParser()
        parser.add_option(
            '-l', '--lamp', action='store',
            dest='arg_lamp', default='',
            help='Output of LAMP')
        parser.add_option(
            '-d', '--dist', action='store',
            dest='arg_dist', default='',
            help='Output of checkPeak.py')
        parser.add_option(
            '-e', '--exp', action='store',
            dest='arg_exp', default='',
            help='Output of checkExp.py')
        parser.add_option(
            '-o', '--out', action='store',
            dest='arg_out', default='',
            help='Output file')
        parser.add_option(
            '-r', '--rank', action='store', type='int',
            dest='arg_rank', default=RANK_THRESHOLD_DEFAULT,
            help='Rank threshold(a positive value)')
        (opt, args) = parser.parse_args()
        arg = opt.__dict__
        lampFile = arg['arg_lamp']
        distFile = arg['arg_dist']
        expFile = arg['arg_exp']
        outFile = arg['arg_out']
        rank_threshold = arg['arg_rank']

        if '' in (lampFile, distFile, expFile, outFile):
            raise TypeError()
        if rank_threshold != RANK_THRESHOLD_DEFAULT and rank_threshold < 1:
            raise TypeError()
    except:
        print 'Usage: ' + str(sys.argv[0]) + \
              ' --lamp out_lamp.txt --dist out_dist.txt' + \
              ' --exp out_exp.txt --out out_report.txt [--rank rank]'
        sys.exit()

    gene2comb, gene2dist, labelArr, combArr = readFiles(
        lampFile, distFile, expFile, rank_threshold)
    o_file = open(outFile, 'w')
    o_file.write('#gene' + SEP + SEP.join(map(str, combArr)) +
                 SEP + SEP.join(map(str, labelArr)) + '\n')
    rlist = sorted(gene2comb.items(), key=lambda x: x[0])

    for ritem in rlist:
        gene = ritem[0]
        o_file.write(gene)

        for comb in combArr:
            value = '-'
            if gene2comb[gene][comb] != DEFAULT_VALUE:
                value = gene2comb[gene][comb]
            o_file.write(SEP + str(value))
        o_file.write(SEP + str(gene2dist[gene]))
    o_file.close()


if __name__ == '__main__':
    main()
