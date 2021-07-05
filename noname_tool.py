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
