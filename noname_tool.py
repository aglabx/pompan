#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.05.2021
#@author: SWW 2021 team
#@contact: support@noname_tool.com

import argparse


def count_GGG_triplets(sequence):
    '''Count GGG triplets in given sequence.
    
    Parameters:
    sequence (string): DNA sequence

    Returns:
    int: A number of GGG triplets
    '''
    return sequence.count('GGG')

import re
import numpy as np
import pprint
from Bio import Seq, SeqUtils, SeqIO

seqs = []
for seq_record in SeqIO.parse("measles.txt", "fasta"):
    seqs.append(seq_record.seq)
    
# DNA codon table
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
def count_codons(seq, dict_codons):
    seq = str(seq)
    for i in range(0, len(seq)-3, 3):
        if seq[i:i+3] not in dict_codons.keys():
            dict_codons[seq[i:i+3]] = 1
        else:
            dict_codons[seq[i:i+3]] += 1
    return dict_codons

# with zero shift
cdns_frq = {}
for id_sq, seq in enumerate(seqs):
    cdns_frq[id_sq] = {}
    cdns_frq[id_sq] = count_codons(seq, cdns_frq[id_sq])
cdns_frq

# with one shift
cdns_frq1 = {}
for id_sq, seq in enumerate(seqs):
    cdns_frq1[id_sq] = {}
    cdns_frq1[id_sq] = count_codons(seq[1:], cdns_frq1[id_sq])
#cdns_frq1

# with two shift
cdns_frq2 = {}
for id_sq, seq in enumerate(seqs):
    cdns_frq2[id_sq] = {}
    cdns_frq2[id_sq] = count_codons(seq[2:], cdns_frq2[id_sq])
#cdns_frq2

# Counting aminas task 
# translating the seq
def translate(seq, protein):
    trans_seq = ''
    for i in range(0, len(seq)-3, 3):
        trans_seq += protein[seq[i: i+3]]  
    return trans_seq

def count_aminas(seq, list_aminas, aminas_dict):

    for amina in list_aminas:
        aminas_dict[amina] = seq.count(amina)
    return aminas_dict

aminas, aminas1, aminas2 = {}, {}, {}
list_aminas = np.unique(list(protein.values()))

# with zero shift
for id_sq, seq in enumerate(seqs):
    aminas[id_sq] = {}
    trans_seq = translate(seq, protein)
    aminas[id_sq] = count_aminas(trans_seq, list_aminas, aminas[id_sq])

    # with one shift
    aminas1[id_sq] = {}
    trans_seq = translate(seq[1:], protein)
    aminas1[id_sq] = count_aminas(trans_seq, list_aminas, aminas1[id_sq])
    
     # with one shift
    aminas2[id_sq] = {}
    trans_seq = translate(seq[2:], protein)
    aminas2[id_sq] = count_aminas(trans_seq, list_aminas, aminas2[id_sq])
    
def main(input_file, output_file):
  
    with open(input_file) as fh:
        pass
    
    for sequence in sequences:
      results = {}
      results['ggg_counts'] = count_GGG_triplets(sequence)
    
    with open(output_file, 'w') as fw:
        pass
  

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Super noname tool!')
    parser.add_argument('-i','--input', help='Input file', required=True)
    parser.add_argument('-o','--output', help='Output file', required=True)
    args = vars(parser.parse_args())
    
    main(input_file, output_file)
