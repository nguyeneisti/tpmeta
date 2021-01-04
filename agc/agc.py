#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Delphine NGUYEN"
__copyright__ = "EISTI / CY Tech"
__credits__ = ["Delphine NGUYEN"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Delphine NGUYEN"
__email__ = "nguyendelp@eisti.eu"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """Retourne un générateur de séquence.
    Séquence comptée si sa longueur est >= à minseqlen"""
    seq_list = []
    with gzip.open(amplicon_file, "rt") as file:
        seq = ""
        for line in file:
            line_s = line.strip("\n")
            if line_s[:1]=="T" or line_s[:1]=="G" or line_s[:1]=="C" or line_s[:1]=="A":
                seq = seq + line_s
            else:
                seq_list.append(seq.strip("\n"))
                seq = ""
        for elm in seq_list:
            if len(elm)>=minseqlen:
                yield elm

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ Retourne un générateur de séquences uniques
    ayant une occurrence >=mincount et leur occurrence
    Format : yield [sequence, count] """
    gen_seq = read_fasta(amplicon_file, minseqlen)
    dict_seq = {}
    for seq in gen_seq:
        if seq in dict_seq:
            dict_seq[seq] += 1
        else:
            dict_seq[seq] = 1
    dict_seq = dict(sorted(dict_seq.items(), key=lambda item:item[1], reverse = True))
    for seq in dict_seq:
        if dict_seq[seq]>=mincount:
            yield [seq, dict_seq[seq]]

def get_chunks(sequence, chunk_size):
    pass

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    pass

def get_identity(alignment_list):
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Retourne un générateur des séquences non chimérique
    Format : yield [sequence, count]"""

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()
