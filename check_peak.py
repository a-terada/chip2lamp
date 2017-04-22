#!/usr/bin/env python
from __future__ import print_function

"""checkPeak.py converts peak bed files to LAMP inputs.

@author:     LAMP dev team
@copyright:  LAMP dev team
@license:    GPLv3
@contact:    lamp_staff(AT)googlegroups.com
@deffield    updated: Updated
"""

import sys
import re
from optparse import OptionParser

__version__ = 1.0
__date__ = '2015-06-27'
__updated__ = '2017-04-22'

# default value
UPDIST_DEFAULT = 2000
INDIST_DEFAULT = 300
SCO_THRESHOLD = 0.0
#DEFAULT_VALUE = -sys.maxint for Python2
#DEFAULT_VALUE = -sys.maxsize for Python3
DEFAULT_VALUE = -1000000


def binFromRangeStandard(sta, end):
    """Calculate bin.

    Keyword arguments:
    sta -- Start position
    end -- End position
    Returns: List, List

    """

    binOffsets = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0]
    binFirstShift = 17
    binNextShift = 3
    startBin = int(sta)
    endBin = int(end) - 1
    startBin >>= binFirstShift
    endBin >>= binFirstShift
    upper = []
    lower = []
    for i in range(len(binOffsets)):
        if startBin == endBin:
            upper.append(binOffsets[i] + startBin)
        else:
            vl = binOffsets[i] + startBin
            vl_2 = binOffsets[i] + endBin
            while vl <= vl_2:
                lower.append(vl)
                vl += 1
        startBin >>= binNextShift
        endBin >>= binNextShift
    return(upper, lower)


def readGeneFile(geneFile, peakFiles):
    """Set default value to gene2peaks.

    Keyword arguments:
    geneFile  -- Gene file in gtf/gff3 format
    peakFiles -- List of peak files
    Returns: Dictionary, Dictionary

    """

    bin2genes = {}
    gene2peaks = {}
    r_geneFile = open(geneFile)
    count = 0
    ecount = 0
    file_part = re.compile('.*\.gtf$', re.I)
    file_part2 = re.compile('.*\.gff3?$', re.I)
    id_part = re.compile('gene_id\s"(\S+)";')
    id_part2 = re.compile('.+Name=(\w+);')
    retsu_part = re.compile('^[-+]?\d+(\.\d+)?([eE][-+]?[0-9]+)?$')

    if (file_part.search(geneFile) is None) \
       and (file_part2.search(geneFile) is None):
        print('Error: Fail to open ' + \
              geneFile + 'with unknown file extentions.',
              file=sys.stderr)
        sys.exit()

    for fh in r_geneFile:
        count += 1
        if fh in ('\n', '\r'):
            ecount += 1
            continue
        if 0 == fh.find('#'):
            if file_part2.search(geneFile) is not None:
                ecount += 1
                continue

        arr = fh.split('\t')
        if 9 > len(arr):
            print('Error: Less columns at line ' + \
                  str(count) + ' in ' + geneFile,
                  file=sys.stderr)
            sys.exit()

        if None is retsu_part.search(arr[3]):
            print(sys.stderr, 'Error: Non-numeric value at line ' + \
                  str(count) + ' column 4 in ' + geneFile, file=sys.stderr)
            sys.exit()

        if None is retsu_part.search(arr[4]):
            print('Error: Non-numeric value at line ' + \
                  str(count) + ' column 5 in ' + geneFile,
                  file=sys.stderr)
            sys.exit()

        chrom = arr[0]
        feat  = arr[2]
        sta = int(arr[3]) - 1
        end = int(arr[4])
        strand = arr[6]
        info = arr[8]
        upperBins, lowerBins = binFromRangeStandard(sta, end)
        ownBin = upperBins[0]
        gene = ''

        if file_part.search(geneFile) is not None:
            if id_part.search(info) is not None and feat == 'exon':
                gene = id_part.match(info).group(1)
        if file_part2.search(geneFile) is not None:
            if id_part2.search(info) is not None \
               and (feat == 'gene' or feat == 'exon'):
                gene = id_part2.match(info).group(1)

        if gene != '':
            if ownBin in bin2genes:
                tmp = bin2genes[ownBin]
            else:
                tmp = []
            x = {'chrom': chrom, 'sta': sta,
                 'end': end, 'strand': strand, 'gene': gene}
            tmp.append(x)
            bin2genes[ownBin] = tmp
            gene2peaks[gene] = {}
            for peakFile in peakFiles:
                gene2peaks[gene].update(
                    {peakFile: {'flag': 0, 'dist': DEFAULT_VALUE}})
    r_geneFile.close()

    if count - ecount <= 0:
        print(sys.stderr, 'Error: No valid line in ' + geneFile,
              file=sys.stderr)
        sys.exit()

    return bin2genes, gene2peaks


def readPeakFile(peakfile, bin2genes, gene2peaks, updist, indist, sco_threshold, macs_flg):
    """Read peak files and set flag & distance to gene2peaks.

    Keyword arguments:
    peakfile -- Peak file
    bin2genes -- Dictionary of bin to genes
    gene2peaks -- Dictionary of gene to peaks
    updist -- Distance upstream from tss
    indist -- Distance downstream from tss
    sco_threshold -- Score threshold
    macs_flg -- True if the peak files are generated by MACS2
    Returns: Dictionary

    """
    col_length = 5
    sco_col = 4 # the columun number for peak score
    if macs_flg:
        col_length = 4
        sco_col = 3
    r_file = open(peakfile, 'r')
    tmp = {}
    count = 0
    ecount = 0
    retsu_part = re.compile('^[-+]?\d+(\.\d+)?([eE][-+]?[0-9]+)?$')
    for fh in r_file:
        count += 1

        if fh in ('\n', '\r'):
            ecount += 1
            continue
        arr = fh.split('\t')
        if col_length > len(arr):
#        if 5 > len(arr):
            print('Error: Less columns at line ' + \
                  str(count) + ' in ' + peakfile, file=sys.stderr)
            sys.exit()

        if None is retsu_part.search(arr[1]):
            print(sys.stderr, 'Error: Non-numeric value at line ' + \
                  str(count) + ' column 2 in ' + peakfile, file=sys.stderr)
            sys.exit()

        if None is retsu_part.search(arr[2]):
            print('Error: Non-numeric value at line ' + \
                  str(count) + ' column 3 in ' + peakfile,
                  file=sys.stderr)
            sys.exit()
        if None is retsu_part.search(arr[sco_col]):
#        if None is retsu_part.search(arr[4]):
            print('Error: Non-numeric value at line ' + \
                  str(count) + ' column 5 in ' + peakfile,
                  file=sys.stderr)
            sys.exit()

        chrom = arr[0]
        sta = int(arr[1])
        end = int(arr[2])
#        sco = float(arr[4])
        sco = float(arr[sco_col])

        if sco < sco_threshold:
            continue
        upperBins, lowerBins = binFromRangeStandard(sta, end)
        nlist = upperBins + lowerBins

        for val in nlist:
            if val not in bin2genes:
                continue
            tmp = bin2genes[val]
            for x in tmp:

                if x['chrom'] != chrom:
                    continue

                prev_dist = gene2peaks[x['gene']][peakfile]['dist']

                # when strand of gene is plus
                if('+' == x['strand']) and \
                        (sta <= (x['sta'] + indist)) and \
                        ((x['sta'] - updist) <= end):
                    gene2peaks[x['gene']][peakfile]['flag'] = 1

                    # 20150303
                    dist = 0
                    if sta < x['sta']:
                      dist = x['sta'] - end
                    else:
                      dist = x['sta'] - sta

                    if(prev_dist == DEFAULT_VALUE) or \
                      (prev_dist > dist):
                        gene2peaks[x['gene']][peakfile]['dist'] = dist

                # strand of gene is minus
                elif('-' == x['strand']) and \
                        (sta <= (x['end'] + updist)) and \
                        ((x['end'] - indist) <= end):
                    gene2peaks[x['gene']][peakfile]['flag'] = 1

                    # 20150303
                    dist = 0
                    if (x['end'] < end):
                      dist = sta - x['end']
                    else:
                      dist = end - x['end']

                    if(DEFAULT_VALUE == prev_dist) or \
                      (prev_dist > dist):
                        gene2peaks[x['gene']][peakfile]['dist'] = dist
    r_file.close()
    if count - ecount <= 0:
        print('Error: No valid line in ' + peakfile, file=sys.stderr)
        sys.exit()

    return gene2peaks


def main():
    try:
        sco_threshold = SCO_THRESHOLD
        peakFiles = []
        geneFile = ''
        distFile = ''
        outFile = ''
        labelStr = ''
        flg = 0
        op_part = re.compile('^-')

        for p in range(len(sys.argv)):
            if sys.argv[p] in ('-p', '--peak'):
                s = p

        for ind, op in enumerate(sys.argv[s + 1:], s + 1):
            if flg == 1:
                continue

            if op_part.search(op) is not None:
                e = ind - 1
                flg = 1
            if ind == len(sys.argv) - 1:
                e = ind
                flg = 1

        n = e - s
        parser = OptionParser()
        parser.add_option(
            '-p', '--peak', action='store', dest='arg_peak',
            default='', nargs=n,
            help='Peak files listed in spaces')
        parser.add_option(
            '-g', '--gene', action='store', dest='arg_gene',
            default='',
            help='Gene file in gtf/gff3 format')
        parser.add_option(
            '-d', '--dist', action='store', dest='arg_dist',
            default='',
            help='Output file (distance)')
        parser.add_option(
            '-o', '--out', action='store', dest='arg_out',
            default='',
            help='Output file (existance)')
        parser.add_option(
            '-u', '--up', action='store', dest='arg_up', type='int',
            default=UPDIST_DEFAULT,
            help='Distance upstream from tss(a non-negative value)')
        parser.add_option(
            '-i', '--in', action='store', dest='arg_in', type='int',
            default=INDIST_DEFAULT,
            help='Distance downstream from tss(a non-negative value)')
        parser.add_option(
            '-l', '--label', action='store', dest='arg_label',
            default='',
            help='Peak names listed in commas')

        # add AT, Apr. 13, 2015
        parser.add_option(
            '--macs2', action='store_true',
            dest='macs2',
            default=False,
            help='Use this option when the peak files are generated with MACS2.')

        (opt, args) = parser.parse_args()
        arg = opt.__dict__
        peakFiles = arg['arg_peak']
        geneFile = arg['arg_gene']
        distFile = arg['arg_dist']
        outFile = arg['arg_out']
        updist = arg['arg_up']
        indist = arg['arg_in']
        labelStr = arg['arg_label']
        if '' in (geneFile, outFile, distFile):
            raise TypeError()

        if 0 == len(peakFiles):
            raise TypeError()

        if updist <= 0 or indist <= 0:
            raise TypeError()

        if str is type(peakFiles):
            peakFiles = [peakFiles]
        else:
            peakFiles = list(peakFiles)
    except:
        print('Usage: ' + str(sys.argv[0]) + \
              ' --gene genes.gtf --peak peakFile [peakFile2 peakFile3 ...] --out out_peak.txt --dist out_dist.txt [--up ' + \
              str(UPDIST_DEFAULT) + '] [--in ' + \
              str(INDIST_DEFAULT) + ']  [--label TF1,TF2, ...] [--macs2]')
        sys.exit()
    check_peak(genefile, peakfiles, tmpoutpeak,
        outdist, updist, indist, labelStr, sco_threshold, arg['macs2'])

def check_peak(geneFile, peakFiles, outFile, distFile,
    updist, indist, labelStr, sco_threshold, macs2):
    if '' == labelStr:
        peakLabels = peakFiles
    else:
        peakLabels = labelStr.split(',')
        # If number of --label is not equal to the number of --peak, discard
        # --label.
        if len(peakLabels) != len(peakFiles):
            peakLabels = peakFiles
    bin2genes, gene2peaks = readGeneFile(geneFile, peakFiles)
    for peakFile in peakFiles:
        gene2peaks = readPeakFile(
            peakFile, bin2genes, gene2peaks, updist, indist,
            sco_threshold, macs2)
    fout = open(outFile, 'w')
    fdist = open(distFile, 'w')

    fout.write('#gene')
    fdist.write('#gene')

    for peakLabel in peakLabels:
        fout.write(',' + peakLabel)
        fdist.write(',' + peakLabel)

    fout.write('\n')
    fdist.write('\n')
    plist = sorted(gene2peaks.items(), key=lambda x: x[0])
    for wline in range(len(plist)):
        fout.write(plist[wline][0])
        fdist.write(plist[wline][0])

        for peakFile in peakFiles:
            flag = str(plist[wline][1][peakFile]['flag'])
            fout.write(',' + flag)
            dist = str(plist[wline][1][peakFile]['dist'])
            if str(DEFAULT_VALUE) in dist:
                dist = '-'
            fdist.write(',' + dist)
        fout.write('\n')
        fdist.write('\n')

    fout.close()
    fdist.close()


if __name__ == '__main__':
    main()
