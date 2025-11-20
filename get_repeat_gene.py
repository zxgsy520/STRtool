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
logging.basicConfig(level=logging.INFO) #初始化时，如果没指定level，那么level的默认级别为WARNING

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    """读取表格文件"""
    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line or line.startswith("#"): #跳过空行和首字母为#的行
            continue

        yield line.split(sep)

    fh.close()


def split_attr(attributes, sep="="):

    """分割gff文件中的最后一列"""
    r = collections.OrderedDict()
    if sep not in attributes:
        sep = " "
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


def read_trf_gff(file):

    r = {}

    for line in read_tsv(file, sep="\t"):
        if line[0] not in r:
            r[line[0]] = []
        attr = split_attr(line[-1])
        r[line[0]].append([int(line[3]), int(line[4]), attr["ID"]])

    return r


def to_string(attributes):

    """输出gff文件的最后一列"""
    attr = []

    for key, value in attributes.items():
        if not value:
            continue
        if key in "ID":
            attr.insert(0, "%s=%s" % (key, value))
        elif key in "Parent":
            attr.insert(1, "%s=%s" % (key, value))
        elif key in "Name":
            attr.insert(1, "%s=%s" % (key, value))
        else:
            attr.append("%s=%s" % (key, value))

    return ";".join(attr)


def get_repeat_gene(gene_gff, trf_gff):
 
    r = read_trf_gff(trf_gff)

    for line in read_tsv(gene_gff, sep="\t"):
        if line[0] not in r:
            continue
        start = int(line[3])
        end = int(line[4])
        temp = r[line[0]]
        for rstart, rend, rid in temp:
            if end < rstart:
                continue
            if start > rend:
                break
            attr = split_attr(line[-1])
            if start >= rstart and end <=end:
                rtype = "include"
            else:
                rtype = "segmentation"
              
            attr["Repeat_id"] = "%s:%s-%s" % (rid, rstart, rend)
            attr["Repeat_type"] = rtype
            line[-1] = to_string(attr)
            print("\t".join(line))

    return 0
       

def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input gene prediction results(genome.gff3)")
    parser.add_argument("--trf_gff", metavar="STR", type=str, default="",
        help="Input TRF prediction results(TRF.gff3)")

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
    get_repeat_gene.py: Obtain genes containing repeat sequences
attention:
    get_repeat_gene.py genome.gff3 --trf_gff genome.TRF.gff3 >gene_TRF.gff3
version: %s
contact:  %s <%s>\
        """ % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    get_repeat_gene(gene_gff=args.input, trf_gff=args.trf_gff)


if __name__ == "__main__":

    main()
