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
            #si la ligne commence par T, G, C ou A, on sauvegarde
            if line_s[:1]=="T" or line_s[:1]=="G" or line_s[:1]=="C" or line_s[:1]=="A":
                seq = seq + line_s
            #autrement on sauvegarde la séquence en entier
            else:
                seq_list.append(seq.strip("\n"))
                seq = ""
        # on renvoie un générateur de séquence
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
    """ Retourne une liste de sous-séquences de taille chunk_size
    non chevauchantes
    Format : [] """
    subseq = []
    for i in range(0, len(sequence), chunk_size):
        if i+chunk_size<=len(sequence):
            subseq.append(sequence[i:i+chunk_size])
    return subseq

def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """Retourne un générateur de tous les mots de longueur k dans sequence
    Format : yield kmer """
    for i in range(0,len(sequence)):
        if i+kmer_size<=len(sequence):
            kmer = sequence[i:i+kmer_size]
            yield kmer
        else:
            break

def get_identity(alignment_list):
    """ Retourne le pourcentage d'identité entre deux séquences """
    count = 0
    length = 0
    length = len(alignment_list[0])
    for j in range(length):
        if alignment_list[0][j] == alignment_list[1][j]:
            count+=1
    return count*100/length

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Retourne générateur séquences non chimérique, Format: yield [seq, count]"""
    gen_seq = []
    for elm in dereplication_fulllength(amplicon_file,minseqlen,mincount):
        gen_seq.append(elm)
    seq_candidate1 = gen_seq[0]
    seq_candidate2 = gen_seq[1]
    seq_candidate1_chunks = get_chunks(seq_candidate1[0], chunk_size)
    seq_candidate2_chunks = get_chunks(seq_candidate2[0], chunk_size)
    dict_simil = {}
    subseq_list = []
    for i in range(1,len(gen_seq)):
        subseq_list.append(get_chunks(gen_seq[i][0], chunk_size))
    for subseq in subseq_list:
        dict_subseq = cut_kmer(subseq, kmer_size)
        for chunk in seq_candidate1_chunks:
            dict_seq1_chunk_candidate = cut_kmer(chunk, kmer_size)
            for chunk2 in seq_candidate2_chunks:
                dict_seq2_chunk_candidate = cut_kmer(chunk2, kmer_size)
                for key in dict_subseq:
                    while len(dict_simil<=8):
                        if key in dict_seq1_chunk_candidate or key in dict_seq2_chunk_candidate:
                            dict_simil[subseq] +=1
    seq_parents = []
    for key in dict_simil:
        if len(seq_parents<=2):
            if key in seq_candidate1 and key in seq_candidate2:
                seq_parents.append(key)
        else:
            break

    identity_list = []
    for elm in seq_parents:
        for i in range(len(get_chunks(elm, chunk_size))):
            alignment_list = nw.global_align(get_chunks(elm[i], chunk_size), seq_candidate1_chunks[i], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identityc1 = get_identity(alignment_list)
            alignment_list = nw.global_align(get_chunks(elm[i], chunk_size), seq_candidate2_chunks[i], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identityc2 = get_identity(alignment_list)
            identity_list.append([identityc1,identityc2])
    for elm in dereplication_fulllength(amplicon_file,minseqlen,mincount):
        yield [elm[0], elm[1]]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ Retourne la liste des OTUS"""
    list = []
    list_identity = []
    chimera = []
    for elm in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
        chimera.append(elm)
    for i in range(len(chimera)):
        alignement_list = nw.global_align(chimera[i][0], chimera[0][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
        identity = get_identity(alignement_list)
        if identity>=0.95:
            list.append([chimera[i][0],chimera[i][1]])
    return list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """ Ecrit dans un fichier externe les OTU"""
    file = open(output_file, "w+")
    for i in range(len(OTU_list)):
        fasta = fill(OTU_list[i][0])
        file.write(">OTU_%d occurrence:%d\n" %(i+1, OTU_list[i][1]))
        file.write(fasta)
        file.write("\n")
    file.close()
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
