#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 21:27:06 2020

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
def calcoor(spec, chrom, start, end):
    size_fp = '/Users/celeste/Documents/Hongkai_lab/dhs_align/align_scripts/chrom_sizes/len/'+ spec +'_chrom_sizes.json'
    with open(size_fp,'r') as g:
        lendict = json.load(g)
    chrom_len = int(lendict[chrom])
    tran_start = chrom_len - end + 1
    tran_end = chrom_len - start + 1
    return tran_start, tran_end

    
# calculate phastcons score (based on human coordinates)
def cal_phastcons_score(path, chrom, start, end):
    with open('/Users/celeste/Documents/Hongkai_lab/dhs_align/align_scripts/chrom_sizes/posi/hg38_chrom_posi.json','r') as g:
        posidict = json.load(g)
    whole_start = int(posidict[chrom]) + start
    whole_end = int(posidict[chrom]) + end
    f = open(path,'rb+')
    f.seek((whole_start-1)*2)
    base_num = whole_end - whole_start + 1
    i = 1
    score = 0
    list1 = []
    while i <= base_num:
        base_score = struct.unpack('H',f.read(2))[0]
        list1.append(base_score)
        score = score + base_score
        i = i + 1
    f.close()
    score = score/(end-start+1)
    return score


def seqdict(ref_seq,cmp_seq):
    seq1 = list(ref_seq.strip().upper())
    seq2 = list(cmp_seq.strip().upper())
    seq_dict = defaultdict()
    for i in range(len(seq1)):
        seq_dict[i] = [seq1[i],seq2[i]]
    return seq_dict


def calpercent(start, end, seqdict):
    same_num = 0
    my_ref_seq = []
    my_cmp_seq = []
    seq_len = 0
    for key, value in seqdict.items():
        if key >= start and key <= end:
            my_ref_seq.append(value[0])
            my_cmp_seq.append(value[1])
            if value[0] == value[1]:
                same_num = same_num + 1
            if value[0] != '-':
                seq_len = seq_len + 1
    percentage = format(same_num/seq_len*100,'0.3f')
    my_ref = ''.join(my_ref_seq)
    my_cmp = ''.join(my_cmp_seq)
    return my_ref, my_cmp, percentage

       
def givcon(midlist,mylist):
    midstart = midlist[0]
    midend = midlist[1]
    mystart = mylist[0]
    myend = mylist[1]
    if int(myend) <= int(midstart):
        return 'low'
    if int(mystart) >= int(midend):
        return 'high'
    if int(myend) > int(midstart) and int(myend) < int(midend) and int(mystart) < int(midstart):
        return 'ovlow'
    if int(mystart) > int(midstart) and int(mystart) < int(midend) and int(myend) > int(midend):
        return 'ovhigh'
    if int(mystart) >= int(midstart) and int(myend) <= int(midend):
        return 'in'
    if int(mystart) <= int(midstart) and int(myend) >= int(midend):
        return 'out'
    

def binarySearch(fp, low, high, region, blockdict):
    if high >= low:
        mid = int(low + (high - low)/2)
        mid_list = blockdict[mid][0].split()[2:4]
        if givcon(mid_list,region) == 'low':
            return(binarySearch(fp, low, mid-1, region,blockdict))
        if givcon(mid_list,region) == 'high':
            return(binarySearch(fp, mid+1, high, region, blockdict))
        if givcon(mid_list,region) in ['ovlow','ovhigh','in','out']:
            return givcon(mid_list,region),mid
    else:
        return 'no match', -1

def contSearch(fp, con, midblock, region, blockdict):
    result = [midblock]
    if con == 'ovlow':
        i =  midblock - 1
        while i >= 0 and i < len(blockdict):
            linelist = blockdict[i][0].split()[2:4]
            if givcon(linelist,region) in ['ovlow','ovhigh','in','out']:
                result.append(i)
            else:
                break
            i = i - 1
    if con == 'ovhigh':
        j = midblock + 1
        while j >= 0 and j < len(blockdict):
            linelist = blockdict[j][0].split()[2:4]
            if givcon(linelist,region) in ['ovlow','ovhigh','in','out']:
                result.append(j)
            else:
                break
            j = j + 1
    if con == 'out':
        a = midblock - 1
        b = midblock + 1
        while a >= 0 and a < len(blockdict) :
            lista = blockdict[a][0].split()[2:4]
            if givcon(lista,region) in ['ovlow','ovhigh','in','out']:
                result.append(a)
            else:
                break
            a = a - 1
        while b >= 0 and b < len(blockdict):
            listb = blockdict[b][0].split()[2:4]
            if givcon(listb,region) in ['ovlow','ovhigh','in','out']:
                result.append(b)
            else:
                break
            b = b + 1
    return result



def outputline(fp, phast, posilist, region, blockdict, rev_spec):
    id_list = []
    myposi_list = []
    outposi_list = []
    percen_list = []
    phast_list = []
    direct_list = []
    for i in posilist:
        theline = blockdict[i][0].split()
        align_id, ref_chrom, ref_start, ref_end, cmp_chrom, cmp_start, cmp_end, align_dire = theline[0],theline[1],int(theline[2]),int(theline[3]),theline[4],int(theline[5]),int(theline[6]),theline[7]
        seq_dict = seqdict(blockdict[i][1],blockdict[i][2])
        my_start, my_end = int(region[0]), int(region[1])
        
        my_posi_ref = ref_start - 1 
        my_posi_cmp = cmp_start - 1
        my_key_start = 0
        my_key_end = len(seq_dict) - 1
        phast_start, phast_end, out_posi_start, out_posi_end = ref_start,ref_end, cmp_start, cmp_end
        for key,value in seq_dict.items():
            if value[0] != '-':
                my_posi_ref = my_posi_ref + 1
            if value[1] != '-':
                my_posi_cmp = my_posi_cmp + 1
            if my_posi_ref == my_start:
                my_key_start = key
                phast_start = my_posi_ref
                out_posi_start = my_posi_cmp
                if value[1] == '-':
                    out_posi_start = my_posi_cmp + 1
            if my_posi_ref == my_end:
                my_key_end = key
                phast_end = my_posi_ref
                out_posi_end = my_posi_cmp
                
        if out_posi_start < out_posi_end: 
            if align_dire == '-':
                out_posi_start, out_posi_end = calcoor(rev_spec, cmp_chrom, out_posi_start, out_posi_end)
                align_dire = '+'
            my_posi = ref_chrom + ':' + str(phast_start) + '-' + str(phast_end)
            out_posi = cmp_chrom + ':' + str(out_posi_start) + '-' + str(out_posi_end)
            out_refseq, out_cmpseq, out_percent = calpercent(my_key_start, my_key_end, seq_dict)
            out_phast = cal_phastcons_score(phast, ref_chrom, phast_start, phast_end)
            
            id_list.append(align_id)
            myposi_list.append(my_posi)
            outposi_list.append(out_posi)
            percen_list.append(out_percent)
            phast_list.append(out_phast)
            direct_list.append(align_dire)
            
    if id_list == []:
        output_line = 'no_match\n'
    else:
        output_line = ''
        for k in range(0,len(id_list)):
            output_line = output_line + myposi_list[k] + '\t' + outposi_list[k]+'\t'+ str(id_list[k])+'\t'+\
                            str(direct_list[k])+'\t'+str(percen_list[k]) +'\t'+ str(phast_list[k])  + '\n'
    return output_line






                
                
    