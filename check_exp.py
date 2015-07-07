#!/usr/bin/env python
"""checkExp.py converts expression files to LAMP input.

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
__updated__ = '2015-06-27'

# default value
Q_THRESHOLD_DEFAULT = 0.05
USE_TYPE_DEFAULT = 'b'
EXP_THRESHOLD_DEFAULT = 0.0


def read_gene_diff_file(gene_file, gene_diff_file,
                        q_threshold, exp_threshold,
                        use_type):
    """generate expression file for LAMP.

    Keyword arguments:
    geneFile -- Gene file in gtf/gff3 format
    geneDiffFile -- Gene expression file created by cuffdiff
    q_threshold -- Maximum threshold of q-value
    exp_threshold -- Minimum threshold of expression value
    use_type -- u: up only, d: down only, b: both
    Returns: Dictionary

    """
    ecount = 0
    count = 0
    gene2exp = {}
    qvals = []
    try:
        fh = open(gene_file, 'r')
    except IOError as e:
        sys.stderr.write(e.strerror + ":" + gene_file + "\n")
        sys.exit()

    file_part = re.compile('.*\.gtf$', re.I)
    file_part2 = re.compile('.*\.gff3?$', re.I)
    id_part = re.compile('gene_id\s"(\S+)";')
    id_part2 = re.compile('.+Name=(\w+);')
    column_part = re.compile('^[-+]?\d+(\.\d+)?([eE][-+]?[0-9]+)?$')
    if (file_part.search(gene_file) is None) \
       and (file_part2.search(gene_file) is None):
        print >> sys.stderr, 'Error: Fail to open ' + \
            gene_file + 'with unknown file extentions.'
        sys.exit()

    for rline in fh:
        count += 1
        if rline in ('\n', '\r'):
            ecount += 1
            continue

        if 0 == rline.find('#'):
            if file_part2.search(gene_file) is not None:
                ecount += 1
                continue

        garr = rline.split('\t')

        if 9 > len(garr):
            print >> sys.stderr, 'Error: Less columns at line ' + \
                str(count) + ' in ' + gene_file
            sys.exit()

        if None is column_part.search(garr[3]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 4 in ' + gene_file
            sys.exit()

        if None is column_part.search(garr[4]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 5 in ' + gene_file
            sys.exit()

        feat = garr[2]
        #if feat == "-":
        #    feat = garr[1]
        info = garr[8]

        if file_part.search(gene_file) is not None:
            if id_part.search(info) is not None and feat == 'exon':
                gene = id_part.match(info).group(1)
                gene2exp[gene] = 0

        if file_part2.search(gene_file) is not None:
            if id_part2.search(info) is not None and \
               (feat == 'gene' or feat == 'exon'):
                gene = id_part2.match(info).group(1)
                gene2exp[gene] = 0

    fh.close()
    if count - ecount <= 0:
        print >> sys.stderr, 'Error: No valid line in ' + gene_file
        sys.exit()

    ecount = 0
    count = 0

    try:
        fh_2 = open(gene_diff_file, 'r')
    except IOError as e:
        sys.stderr.write(e.strerror + ":" + gene_diff_file + "\n")
        sys.exit()

    for dline in fh_2:
        count += 1
        if (0 == dline.find('test_id')) or (dline in ('\n', '\r')) or \
           (0 == dline.find('#')):
            ecount += 1
            continue
        darr = dline.split('\t')
        if 13 > len(darr):
            print >> sys.stderr, 'Error: Less columns at line ' + \
                str(count) + ' in ' + gene_diff_file
            sys.exit()

        if None is column_part.search(darr[7]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 8 in ' + gene_diff_file
            sys.exit()

        if None is column_part.search(darr[8]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 9 in ' + gene_diff_file
            sys.exit()

        if None is column_part.search(darr[12]):
            print >> sys.stderr, 'Error: Non-numeric value at line ' + \
                str(count) + ' column 13 in ' + gene_diff_file
            sys.exit()

        gene_str = darr[2]
        gene_arr = gene_str.split(',')
        val1 = float(darr[7])
        val2 = float(darr[8])
        q = float(darr[12])
        #sys.stderr.write("%s %f\n" % (gene_str, q))

        if val1 < val2:
            tipe = 'u'
        else:
            tipe = 'd'

        if gene_str == '-':
            continue
        if q > q_threshold:
            continue
        if (use_type != 'b') and (use_type != tipe):
            continue
        if (val1 < exp_threshold) and (val2 < exp_threshold):
            continue

        for gline in gene_arr:
            if gline in gene2exp:
                gene2exp[gline] = 1
                qvals.append(q)
    fh_2.close()
    if count - ecount <= 0:
        print >> sys.stderr, 'Error: No valid line in ' + gene_diff_file
        sys.exit()
    rlist = sorted(gene2exp.items(), key=lambda x: x[0])
    return rlist

def check_exp(gene_file, gene_diff_file,
               q_threshold, exp_threshold,
               use_type, out_file):
    gene2exp = read_gene_diff_file(gene_file, gene_diff_file,
                   q_threshold, exp_threshold,
                   use_type)
    fout = open(out_file, 'w')
    fout.write('#gene,expression' + '\n')
    for wline in range(len(gene2exp)):
        fout.write(str(gene2exp[wline][0]) + ',' +
                   str(gene2exp[wline][1]) + '\n')
    fout.close()


def main():
    gene2exp = {}
    use_type = USE_TYPE_DEFAULT

    try:
        parser = OptionParser()
        parser.add_option(
            '-g', '--gene', action='store', dest='arg_gene',
            default='',
            help='Gene file in gtf/gff3 format')
        parser.add_option(
            '-d', '--diff', action='store', dest='arg_diff',
            default='',
            help='Gene expression file created by cuffdiff')
        parser.add_option(
            '-o', '--out', action='store', dest='arg_out',
            default='',
            help='Output file')
        parser.add_option(
            '-q', '--qval', action='store', dest='arg_qval', type='float',
            default=Q_THRESHOLD_DEFAULT,
            help='Maximum thhreshold of q-value(a non-negative value)')
        parser.add_option(
            '-e', '--exp', action='store', dest='arg_exp', type='float',
            default=EXP_THRESHOLD_DEFAULT,
            help='Minimum threshold of expression value(a non-negative value)')
        (opt, args) = parser.parse_args()
        arg = opt.__dict__
        gene_file = arg['arg_gene']
        gene_diff_file = arg['arg_diff']
        out_file = arg['arg_out']
        q_threshold = arg['arg_qval']
        exp_threshold = arg['arg_exp']
        if '' in (gene_file, gene_diff_file, out_file):
            raise TypeError()
        if q_threshold < 0 or exp_threshold < 0:
            raise TypeError()
    except:
        print 'Usage:' + ' ' + str(sys.argv[0]) + \
            ' --gene genes.gtf --diff gene_exp.diff --out out_exp.txt [--qval ' + \
              str(Q_THRESHOLD_DEFAULT) + '] [--exp ' + \
              str(EXP_THRESHOLD_DEFAULT) + ']'
        sys.exit()

    check_exp(gene_file, gene_diff_file,
               q_threshold, exp_threshold,
               use_type, out_file)


if __name__ == '__main__':
    main()
