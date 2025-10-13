#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.1.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def get_alias(repeat_unit):

    unit_len = len(repeat_unit)
    alias = ""
    if unit_len == 1:
        alias = "monomer"
    elif unit_len == 2:
        alias = "dimer"
    elif unit_len == 3:
        alias = "trimer"
    elif unit_len == 4:
        alias = "tetramer"
    elif unit_len == 5:
        alias = "pentamer"
    elif unit_len == 6:
        alias = "hexamer"
    elif unit_len == 7:
        alias = "heptamer"
    elif unit_len == 8:
        alias = "octamer"
    elif unit_len == 9:
        alias = "nonamer"
    elif unit_len == 10:
        alias = "decamer"
    else:
        alias = "%smer" % unit_len

    return alias


def read_trf_dat(file, repeat_num=2, percent_matches=100.0, sep=None):

    if repeat_num:
        repeat_num = float(repeat_num)

    r = []

    source = "TRF:4.04"
    seqid = ""
    type = "tandem_repeat"
    phase = "."
    score = "."

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue
        if line.startswith("Version"):
            source = "TRF:%s" % line.split()[-1]
            continue
        if line.startswith("Sequence"):
            seqid = line.split()[1]
            continue
        line = line.split(sep)

        if len(line) <=12:
            continue
        start, end = line[0], line[1]
        strand = "+"
        if int(start) >= int(end):
            strand = "-"
            start, end = end, start
        alias = get_alias(line[13])
        rnum = line[3]
        matche = line[5]
        attributes = "Alias=%s;Name=%s;Repeat_num=%s;Percent_matches=%s" % (alias, line[13], rnum, matche)

        if repeat_num:
            if repeat_num != float(rnum): #过滤重复次数
                continue
        if float(matche) < percent_matches: #过滤相似度低的序列
            continue

        r.append([seqid, source, type, start, end, score, strand, phase, attributes])

    return r


def split_attr(attributes):

    r = collections.OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if "=" not in content:
            LOG.info("%r is not a good formated attribute: no tag!")
            continue
        tag, value = content.split("=", 1)
        r[tag] = value

    return r


def _overlap(q, s):

    s1, e1 = int(q[3]), int(q[4])
    s2, e2 = int(s[3]), int(s[4])
    assert s2 >= s1

    if (e1-s2+1) > min((e1-s1+1, e2-s2+1))*0.6:
        return 1
    else:
        return 0


def deduplicate_trf(gfflist):

    r = []
    p = []

    for g in sorted(gfflist, key=lambda i: (i[0], int(i[3]), int(i[4]))):

        if not p:
            p = g
            continue

        if g[0] != p[0]:
            r.append(p)
            p = g
            continue

        if _overlap(p, g):
            if int(g[4])-int(g[3]) > int(p[4])-int(p[3]):
                p = g
            continue
        r.append(p)
        p = g

    if p:
        r.append(p)

    return r


def trf2gff(files, repeat_num=2, percent_matches=100.0, outfa=""):

    r = []
    for file in files:
        r += read_trf_dat(file, repeat_num, percent_matches)

    if outfa:
        fa = open(outfa, "w")

    n = 0
    data = {}
    total_trf = 0
    total_number = 0
    temp = []

    for line in deduplicate_trf(r):
        n += 1
        attr  = split_attr(line[-1])
        line[-1] = "ID=TRF%s;%s" % (n, line[-1])
        print("\t".join(line)) #输出重复序列的gff

        if outfa:
            fa.write(">TRF%s\n%s\n" % (n, attr["Name"]))

        unit_len = len(attr["Name"])
        if unit_len not in data:
            data[unit_len] = []
        trflen = int(line[4])-int(line[3])+1
        data[unit_len].append(trflen)
        total_trf += trflen
        total_number += 1

    if outfa:
        fa.close()  
    fo = open("stat_trf.tsv", "w")
    fo.write("#Motif(-mer)\tCount\tLength\tCount Percentage(%)\n")
    temp_number = 0
    temp_len = 0
    for i, lengths in sorted(data.items(), key=lambda x: x[0]):
        sum_count = len(lengths)
        temp_number += sum_count
        temp_len += sum(lengths)

        if i >= 16:
            break
        fo.write("{0}\t{1:,}\t{2:,}\t{3:.2f}\n".format(i, sum_count, sum(lengths), sum_count*100.0/total_number))

    fo.write("{0}\t{1:,}\t{2:,}\t{3:.2f}\n".format("Other", total_number-temp_number, 
        total_trf-temp_len, (total_number-temp_number)*100.0/total_number)
    )
    fo.write("{0}\t{1:,}\t{2:,}\t{3:.2f}\n".format("Total", total_number, total_trf, 100))
    fo.close()

    return 0

def add_hlep_args(parser):

    parser.add_argument("input", nargs="+", metavar="FILE", type=str,
        help="Input TRF prediction results, (*.dat)")
    parser.add_argument("-rn", "--repeat_num", metavar="STR", type=str, default="",
        help="Set the filter repeat number, default=")
    parser.add_argument("-pm", "--percent_matches", metavar="FLOAT", type=float, default=50,
        help="Set the percentage of filter matches, default=50") 
    parser.add_argument("-of", "--outfa", metavar="FILE", type=str, default="",
        help="Set output repeat sequence, default=")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""
name:
    trf2gff.py: Merge and convert TRF duplicate prediction result files
attention:
    trf2gff.py *.TRF.dat >TRF.gff3
version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    trf2gff(files=args.input, repeat_num=args.repeat_num, percent_matches=args.percent_matches, outfa=args.outfa)


if __name__ == "__main__":

    main()
