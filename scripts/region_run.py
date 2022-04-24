#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 17:06:55 2021

@author: celeste
"""
import region_modules as rm
from collections import defaultdict

workdir = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/'
inputfile = '/Users/celeste/Documents/Hongkai_lab/dhs_align/mm10_hg38/ENCFF910SRW.tsv'
o = open(workdir + 'mm10_hg38_dhs.txt','a')
error = []
#phast_name = workdir + 'hg38_phastConsElements100way.bin'

# coordiantes transformation if alignment on negative strand
rev_spec = 'hg38'
# axt net file suffix
axt_fp = workdir + 'axtnet/'
axt_suffix = '.mm10.hg38.net.axt'

chrom_0 = 'chr0'
with open(inputfile,'r') as i:
    for line in i:
        try:
            linelist = line.strip().split('\t')
            index = linelist[3]
            chrom = linelist[0]
            start = int(linelist[1])
            end = int(linelist[2]) 
            mylist = [str(start),str(end)]
            refile = axt_fp + chrom + axt_suffix
            
            if chrom != chrom_0:
                block_index = 0
                block_dict =  defaultdict()           
                for each in rm.readblock(refile):
                    block_dict[block_index] =  each[:]
                    block_index = block_index + 1
            
            con,posi = rm.binarySearch(refile, 0, len(block_dict)-1, mylist, block_dict)
            if con == 'in':
                out = rm.outputline(refile, [posi], mylist, block_dict, rev_spec)
                if out != 'no_match\n':
                    final_line = index + '\t' + out
                    o.write(final_line)
            if con in ['ovlow','ovhigh','out']:
                posilist = rm.contSearch(refile, con, posi, mylist, block_dict)
                out = rm.outputline(refile, posilist, mylist, block_dict, rev_spec)
                if out != 'no_match\n':
                    final_line = ''
                    for eachline in out.split('\n'):
                        if eachline != '':
                            final_line = final_line + index + '\t' + eachline + '\n'
                    o.write(final_line)
            chrom_0 = chrom
        except Exception as e:
            print(e)
            print(line)
            error.append(index)
o.close()

    
    
    
    
