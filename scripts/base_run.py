#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 18:00:47 2020

@author: celeste
"""

import base_modules as bm
from collections import defaultdict


workdir = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/'
inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/ENCFF910SRW_summit.tsv'

#phast_name = workdir + 'hg38_phastConsElements100way.bin'
o = open(workdir + 'mm10_hg38_dhs_summit_extended.txt','a')
error = []

# coordiantes transformation if alignment on negative strand
rev_spec = 'hg38'
# axt net file suffix
axt_fp = workdir + 'axtnet/'
axt_suffix = '.mm10.hg38.net.axt'
# extent region length
setlen = 200

chrom_0 = 'chr0'
with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[3]
            chrom = linelist[0]
            start = int(linelist[1])
            end = int(linelist[2])
            base = int(linelist[4])
            refile = axt_fp + chrom + axt_suffix
            
            if chrom != chrom_0:
                block_index = 0
                block_dict =  defaultdict()           
                for each in bm.readblock(refile):
                    block_dict[block_index] =  each[:]
                    block_index = block_index + 1
                    
            posi = bm.binarySearch(refile, 0, len(block_dict)-1, base, block_dict)
            if posi != 'no match':
                out = bm.extend_outputline(refile, posi, base, setlen, block_dict, rev_spec)
                if out != 'no match\n':
                    final_line = index + '\t' + out
                    o.write(final_line)
            chrom_0 = chrom
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()
