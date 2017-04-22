#!/usr/bin/env python
"""chip2lamp.py convert outputs of MACS and Cuffdiff to the input of LAMP.

@author:     LAMP dev team
@copyright:  LAMP dev team
@license:    GPLv3
@contact:    lamp_staff(AT)googlegroups.com
@deffield    updated: Updated

"""

import sys
import os
import re
import subprocess
from optparse import OptionParser
from check_exp import check_exp
from check_peak import check_peak

# __all__ = []
__version__ = 1.0
__date__ = '2015-06-27'
__updated__ = '2017-04-22'

DEBUG = 0
TESTRUN = 0
PROFILE = 0

# default value
UPDIST_DEFAULT = 2000
INDIST_DEFAULT = 300
Q_THRESHOLD_DEFAULT = 0.05
EXP_THRESHOLD_DEFAULT = 0.0
OUT_DEFAULT = 'out'


def usage(program_name):
    """Print usage.

    Keyword arguments:
    program_name -- Program name

    Returns: None

    """
    print('Usage: ' + program_name + ' --gene genes.gtf --diff gene_exp.diff --peak peakFile [peakFile2 peakFile3 ...]' \
        ' [--exp ' + str(EXP_THRESHOLD_DEFAULT) + '] [--qval ' + str(Q_THRESHOLD_DEFAULT) + ']' \
        ' [--up ' + str(UPDIST_DEFAULT) + '] [--in ' + str(INDIST_DEFAULT) + '] [--out ' + str(OUT_DEFAULT) + ']' \
        ' [--label TF1,TF2, ...] [--peakcheck] [--macs2]')


def checkAllZero(arg0):
    """Check if values in list are all zero.

    Keyword arguments:
    arg0 -- Numeric list

    Returns: True or False

    """
    retVal = True
    for dt in arg0:
        if dt != '0':
            retVal = False
    return retVal


def check_consistency(expfile, peakfile, outexpfile, outpeakfile, peak_check):
    """Check consistency between outputs of checkPeak.py and checkExp.py.

    Keyword arguments:
    expfile -- Temporary output file of checkExp.py
    peakfile -- Temporary output file (peak) of checkPeak.py
    outexpfile -- Final output file of checkExp.py
    outpeakfile -- Final output file of checkPeak.py
    peak_check -- 0: Use all genes, 1: Discard genes not binding to any TF.

    Returns: None

    """
    peaks = {}
    sp = []
    headprog = re.compile('^#')

    with open(outpeakfile, 'w') as outpeakfp, open(outexpfile, 'w') as outexpfp:
        num = 0
        for line in open(peakfile, 'r'):
            num = num + 1
            line = line.rstrip()
            if headprog.match(line):
                print(line, file=outpeakfp)
                continue
            sp = line.split(',')
            size = len(sp)
            if size < 2:
                print('Error: Few columns at line %d in %s.' % (
                    num, peakfile), file=sys.stderr)
                sys.exit()
            peaks[sp[0]] = [checkAllZero(sp[1:]), line]

        num = 0
        for line in open(expfile, 'r'):
            num = num + 1
            line = line.rstrip()
            if headprog.match(line):
                print(line, file=outexpfp)
                continue
            sp = line.split(',')
            size = len(sp)
            if size < 2:
                print('Error: Few columns at line %d in %s.' % (
                    num, peakfile), file=sys.stderr)
                sys.exit()
            if peaks.get(sp[0], '') == '':
                print('Warning: Gene %s is not found in %s.' % (
                    sp[0], peakfile), file=sys.stderr)
                continue

            if not peak_check:
                outexpfp.write(str(line) + '\n')
                outpeakfp.write(str(peaks[sp[0]][1]) + '\n')
            elif not peaks[sp[0]][0]:
                outexpfp.write(str(line) + '\n')
                outpeakfp.write(str(peaks[sp[0]][1]) + '\n')
            peaks[sp[0]] = ''

    for (key, value) in peaks.items():
        if value is not '':
            print('Warning: Gene %s is not found in %s.' % (
                key, expfile), file=sys.stderr)

def main(argv=None):
    program_name = os.path.basename(sys.argv[0])

    if argv is None:
        argv = sys.argv[1:]
    #try:
        try:
            flg = 0
            s = 0
            e = 0
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
                '-g', '--gene', action='store',
                dest='arg_gene',
                default='',
                help='Gene file in gtf/gff3 format')
            parser.add_option(
                '-d', '--diff', action='store',
                dest='arg_diff',
                default='',
                help='Gene expression file created by cuffdiff')
            parser.add_option(
                '-o', '--out', action='store',
                dest='arg_out',
                default=OUT_DEFAULT,
                help='Output prefix')
            parser.add_option(
                '-q', '--qval', action='store',
                dest='arg_qval',
                default=Q_THRESHOLD_DEFAULT,
                help='Maximum thhreshold of q-value(a non-negative value)')
            parser.add_option(
                '-x', '--exp', action='store',
                dest='arg_exp',
                default=EXP_THRESHOLD_DEFAULT,
                help='Minimum threshold of expression value(a non-negative value)')
            parser.add_option(
                '-p', '--peak', action='store',
                dest='arg_peak',
                default='', nargs=n,
                help='Peak files listed in spaces')
            parser.add_option(
                '-u', '--up', action='store',
                dest='arg_up',
                default=UPDIST_DEFAULT,
                help='Distance upstream from tss(a non-negative value)')
            parser.add_option(
                '-i', '--down', action='store',
                dest='arg_in',
                default=INDIST_DEFAULT,
                help='Distance downstream from tss(a non-negative value)')
            parser.add_option(
                '-l', '--label', action='store',
                dest='arg_label',
                default='',
                help='Peak names listed in commas')
            parser.add_option(
                '-c', '--nm-ignore', action='store_true',
                dest='arg_peakcheck',
                default=False,
                help='Check peak(True or False)')
            parser.add_option(
                '-v', '--verbose', action='store',
                dest='arg_verbose',
                default=False,
                help='v')

            # add AT, Apr. 13, 2015
            parser.add_option(
                '--macs2', action='store_true',
                dest='macs2',
                default=False,
                help='Use this option when the peak files are generated with MACS2.')

            (o, a) = parser.parse_args()
            op = o.__dict__
            genefile = op['arg_gene']
            difffile = op['arg_diff']
            out = op['arg_out']
            q_threshold_default = float(op['arg_qval'])
            exp_threshold_default = float(op['arg_exp'])
            peakfiles = op['arg_peak']
            updist = int(op['arg_up'])
            indist = int(op['arg_in'])
            peak_check = op['arg_peakcheck']
            verbose = op['arg_verbose']
            label = op['arg_label']
            macs2 = op['macs2']
            if '' in (genefile, difffile):
                raise TypeError()

            if 0 == len(peakfiles):
                raise TypeError()
        except:
            usage(program_name)
            return 2

        outexp  = out + '_exp.txt'
        outpeak = out + '_peak.txt'
        outdist = out + '_dist.txt'

        print("Upstream from TSS (bp): %d" % updist, file=sys.stderr)
        print("Downstream from TSS (bp): %d" % indist, file=sys.stderr)
        print("q-value threshold for DEG: %f" % q_threshold_default, file=sys.stderr)

        # Execute checkExp.pl
        tmpoutexp = outexp + '.tmp'
        check_exp(genefile, difffile,
                q_threshold_default,
                exp_threshold_default,
                'b', tmpoutexp)

        # Execute checkPeak.pl
        tmpoutpeak = outpeak + '.tmp'
        check_peak(genefile, peakfiles, tmpoutpeak,
            outdist, updist, indist, label, 0.0, macs2)

        check_consistency(tmpoutexp, tmpoutpeak, outexp, outpeak, peak_check)

        # Remove temporary files
        if not verbose:
            os.remove(tmpoutexp)
            os.remove(tmpoutpeak)

    #except Exception, e:
    #    indent = len(program_name) * ' '
    #    sys.stderr.write(program_name + ': ' + repr(e) + '\n')
    #    sys.stderr.write(indent + '  for help use --help\n')
    #    return 2


if __name__ == '__main__':
    if DEBUG:
        sys.argv.append('-h')
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'checkConsistency_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open('profile_stats.txt', 'wb')
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
