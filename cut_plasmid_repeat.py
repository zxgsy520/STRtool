#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import gzip
import logging
import argparse

import collections

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)

    fh.close()


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    r = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if r:
                yield r.split("\n", 1)
            r = "%s\n" % line.strip(">") #去除大于号，没有去除末尾的注释
            continue
        r += line.upper()

    if r:
        yield r.split("\n", 1)
    fp.close()


def split_attr(attributes, sep="="):

    """分割gff文件中的最后一列"""
    r = collections.OrderedDict()
    contents = attributes.split(";")

    for content in contents:
        if not content:
            continue
        if sep not in content:
            LOG.info("%r is not a good formated attribute: no tag!" % content)
            continue
        tag, value = content.strip().split(sep, 1)
        r[tag] = value.strip('"')

    return r


def cut_plasmid_repeat(input, gff, prefix="out"):

    r = {}

    for line in read_tsv(gff, "\t"):
        if line[0] not in r:
            r[line[0]] = []
        attr = split_attr(line[-1])
        temp = [int(line[3]), int(line[4]), attr["Name"]]
        r[line[0]].append(temp)

    fo = open("%s.cut_TRF.bed" % prefix, "w")
    fo.write("Seqid\tStart\tEnd\tCut-Start\tCut-End\n")
    for seqid, seq in read_fasta(input):
        seqid = seqid.split()[0]
        if seqid not in r:
            continue
        temp = r[seqid]
        t1 = 1
        for start, end, rqseq in sorted(temp, key=lambda x: x[0]):
            n_start = start-t1
            n_end = start+len(rqseq)-t1 #参考序列没有错，主要是位置坐标可能不对
            n_rqseq = seq[n_start:n_end]
            n_id = "%s-%s_%s" % (seqid, start, start+2*len(rqseq))
            #print(">%s-%s_%s\n%s" % (seqid, start, start+2*len(rqseq), rqseq))
            fo.write("%s\t%s\t%s\t%s\t%s\n" % (seqid, start, start+2*len(rqseq), n_start+1, n_end-1)) 
            if not seq[n_end::].startswith(rqseq):
                #LOG.info("验证成功:%s\t%s" % (n_id, n_rqseq))
                LOG.info("Validation failed: %s\t%s" % (n_id, n_rqseq))
            seq = seq[0:n_start] + seq[n_end::]
            t1 = t1+len(rqseq)
        print(">%s\n%s" % (seqid, seq))

    fo.close()


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input genome sequence(fasta)")
    parser.add_argument("-g", "-gff", metavar="FILE", type=str, required=True,
        help="Input genome TRF prediction results, 77mer-TRF.gff3")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="out",
        help="Output file prefix")

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
    cut_plasmid_repeat.py: Cut the TRF repeat sequence of the plasmid
attention:
    cut_plasmid_repeat.py plasmid.77mer-TRF.fasta --gff plasmid.77mer-TRF.gff3 >plasmid.cut_77mer-TRF.fasta 2>log
version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    cut_plasmid_repeat(input=args.input, gff=args.gff, prefix=args.prefix)


if __name__ == "__main__":

    main()

