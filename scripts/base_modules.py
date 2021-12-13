#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 18:06:38 2020

@author: celeste
"""

import struct
import json
from collections import defaultdict

# generator of each alignment block in axtnet file
def readblock(fp):
    lines = []
    with open(fp) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line != '\n':
                lines.append(line.strip())
            if line == '\n':
                yield lines
                lines.clear()
        if lines!= []:
            yield line


# transform coordinates if on negative strand
def calcoor(spec, chrom, base):
    size_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/align_scripts/chrom_sizes/len/'+ spec +'_chrom_sizes.json'
    with open(size_fp,'r') as g:
        lendict = json.load(g)
    chrom_len = int(lendict[chrom])
    tran_base = chrom_len - base + 1
    return tran_base

    
# calculate phastcons score (based on human coordinates)
# def cal_phastcons_score(path, chrom, start, end):
#     with open('/Users/celeste/Documents/Hongkai_lab/dhs_align/align_scripts/chrom_sizes/posi/hg38_chrom_posi.json','r') as g:
#         posidict = json.load(g)
#     whole_start = int(posidict[chrom]) + start
#     whole_end = int(posidict[chrom]) + end
#     f = open(path,'rb+')
#     f.seek((whole_start-1)*2)
#     base_num = whole_end - whole_start + 1
#     i = 1
#     score = 0
#     list1 = []
#     while i <= base_num:
#         base_score = struct.unpack('H',f.read(2))[0]
#         list1.append(base_score)
#         score = score + base_score
#         i = i + 1
#     f.close()
#     score = score/(end-start+1)
#     return score


def seqdict(ref_seq,cmp_seq):
    seq1 = list(ref_seq.strip().upper())
    seq2 = list(cmp_seq.strip().upper())
    seq_dict = defaultdict()
    for i in range(len(seq1)):
        seq_dict[i] = [seq1[i],seq2[i]]
    return seq_dict

def getseq(start, end, seqdict):
    my_ref_seq = []
    my_cmp_seq = []
    for key, value in seqdict.items():
        if key >= start and key <= end:
            my_ref_seq.append(value[0])
            my_cmp_seq.append(value[1])
    return my_ref_seq, my_cmp_seq

def searchclosest(key, seqdict):
    start = key
    end = key
    while start >= 0:
        if seqdict[start][1] != '-':
            mystart = start
            break
        start = start - 1
    while end <= len(seqdict)-1:
        if seqdict[end][1] != '-':
            myend = end
            break
        end = end + 1
    inter1 = key - mystart
    inter2 = myend - key
    if inter1 <= inter2:
        return 'before', mystart
    if inter1 > inter2:
        return 'after', myend
       
def givcon(midlist,base):
    midstart = int(midlist[0])
    midend = int(midlist[1])
    mybase = int(base)
    if mybase < midstart:
        return 'low'
    if mybase > midend:
        return 'high'
    if mybase >= midstart and mybase <= midend:
        return 'in'

def binarySearch(fp, low, high, base, blockdict):
    if high >= low:
        mid = int(low + (high - low)/2)
        mid_list = blockdict[mid][0].split()[2:4]
        if givcon(mid_list,base) == 'low':
            return(binarySearch(fp, low, mid-1, base, blockdict))
        if givcon(mid_list,base) == 'high':
            return(binarySearch(fp, mid+1, high, base, blockdict))
        if givcon(mid_list,base) == 'in':
            return mid
    else:
        return 'no match'
    

def base_outputline(fp, posi, base, setlen, blockdict, rev_spec):
    theline = blockdict[posi][0].split()
    align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]),theline[7]
    seq_dict = seqdict(blockdict[posi][1],blockdict[posi][2])

    my_base = int(base)
    my_posi_ref = ref_start - 1
    my_posi_cmp = cmp_start - 1
    out_posi = cmp_start
    for key,value in seq_dict.items():
        if value[0] != '-':
            my_posi_ref = my_posi_ref + 1
        if value[1] != '-':
            my_posi_cmp = my_posi_cmp + 1
        if my_posi_ref == my_base:
            out_posi = my_posi_cmp
            my_key = key
            status = 'non-gap'
            if value[1] == '-':
                status = 'gap'
                closest, thekey = searchclosest(key, seq_dict)
                if closest == 'after':
                    out_posi = out_posi + 1
                    my_key = key + 1
            break
    ref_seq, cmp_seq = getseq(my_key-setlen/2+1, my_key+setlen/2, seq_dict)
    if cmp_seq.count('-')/setlen <= 0.5:
        if align_dire == '-':
            out_posi = calcoor(rev_spec, cmp_chrom, out_posi)
            align_dire = '+'
        my_posi = ref_chrom + ':' + str(base) 
        out_posi = cmp_chrom + ':' + str(out_posi)
        output_line = my_posi + '\t' + out_posi + '\t' + str(align_id) + '\t'+ align_dire +'\t'+ status + '\n'
    else:
        output_line = 'no match\n'
    return output_line, cmp_seq.count('-')

def extend_outputline(fp, posi, base, setlen, blockdict, rev_spec):
    theline = blockdict[posi][0].split()
    align_id, ref_chrom, ref_start, cmp_chrom, cmp_start, align_dire = theline[0],theline[1],int(theline[2]),theline[4],int(theline[5]),theline[7]
    seq_dict = seqdict(blockdict[posi][1],blockdict[posi][2])

    my_base = int(base)
    my_posi_ref = ref_start - 1
    my_posi_cmp = cmp_start - 1
    out_posi = cmp_start
    for key,value in seq_dict.items():
        if value[0] != '-':
            my_posi_ref = my_posi_ref + 1
        if value[1] != '-':
            my_posi_cmp = my_posi_cmp + 1
        if my_posi_ref == my_base:
            out_posi = my_posi_cmp
            my_key = key
            status = 'non-gap'
            if value[1] == '-':
                status = 'gap'
                closest, thekey = searchclosest(key, seq_dict)
                if closest == 'after':
                    out_posi = out_posi + 1
                    my_key = key + 1
            break
    ref_seq, cmp_seq = getseq(my_key-int(setlen/2)+1, my_key+int(setlen/2), seq_dict)
    if cmp_seq.count('-')/setlen <= 0.5:
        if align_dire == '-':
            out_posi = calcoor(rev_spec, cmp_chrom, out_posi)
            align_dire = '+'
        my_posi = ref_chrom + ':' + str(base-int(setlen/2)+1) + '-' + str(base+int(setlen/2))
        out_posi = cmp_chrom + ':' + str(out_posi-int(setlen/2)+1) + '-' + str(out_posi+int(setlen/2))
        output_line = my_posi + '\t' + out_posi + '\t' + str(align_id) + '\t'+ align_dire +'\t'+ status + '\n'
    else:
        output_line = 'no match\n'
    return output_line, cmp_seq.count('-')
    
