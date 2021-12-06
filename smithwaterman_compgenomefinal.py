#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:00:35 2021

@author: eleanorhilgart
"""
import sys
sys.path.append("/Users/eleanorhilgart/github/cg-66-contaminant-detection/data_utils")


import numpy
import time

from data_utils import read_fasta_files
from data_utils import read_fastq_files
from fasta import FASTA
from fastq import FASTQ

def cost(xc, yc):
# Cost function
    if xc == yc: return 4 # match
    if xc == '-' or yc == '-': return -6 # gap
    return -4

def smithWaterman(x, y, s):
# Smith Waterman from Lecture
    V = numpy.zeros((len(x)+1, len(y)+1), dtype=int)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            V[i, j] = max(V[i-1, j-1] + s(x[i-1], y[j-1]), 
                          V[i-1, j  ] + s(x[i-1], '-'),    
                          V[i  , j-1] + s('-',    y[j-1]), 
                          0)                             
    maxlastrow = V.max(axis = 1)[-1] 
    #changed to get the max of the last row as required for read alignment
    return V, maxlastrow

def sw_overfile(readfh, humfhs, contamfhs):
    ''' Iterate over reads, comparing to references (provided in lists), and return dictionary 
    of read ids sorted by whether the reads map to human, contaminant, or ambiguous '''
    # read reference files
    humFASTAs = read_fasta_files(humfhs)
    contamFASTAs = read_fasta_files(contamfhs)
    fq = read_fastq_files(readfh)


    detection = {'Contaminant': [], 'Human': [], 'Cannot be assigned': []}
    #temporary until we decide on standard output format for detection
    
    for val in fq:
        seqid = val['id']
        
        contamvals = []
        humvals = []

        for conFA in contamFASTAs:
            contam = [conFA.get_data()]
            for a in contam:
                contamSW, topcontamval = smithWaterman(val['seq'], a, cost)
                contamvals.append(topcontamval)
        
        for humFA in humFASTAs:
            chrom = [humFA.get_data()]
            for b in chrom:
                humSW, tophumval = smithWaterman(val['seq'], b, cost)
                humvals.append(tophumval)
            
        bestcontam = max(contamvals)
        besthum = max(humvals)
        
        if bestcontam > besthum:
            detection['Contaminant'].append(seqid)
        elif bestcontam < besthum:
            detection['Human'].append(seqid)
        elif bestcontam == besthum:
            detection['Cannot be assigned'].append(seqid)

    return detection

start_time = time.time()


reads = '/Users/eleanorhilgart/Desktop/compgenom_data/compgenom_cleanedfinaldata/simulated_reads/small:tiny data/chr1_Mfermentans_8020_hiseq_reads_tinycut_R1.fastq'
hum = ['/Users/eleanorhilgart/Desktop/compgenom_data/compgenom_cleanedfinaldata/compgenome_readsim_fastas/noN_chr1_cut.fasta']
contam = ['/Users/eleanorhilgart/Desktop/compgenom_data/compgenom_cleanedfinaldata/compgenome_readsim_fastas/Mfermentansbac1_cut.fasta', '/Users/eleanorhilgart/Desktop/compgenom_data/compgenom_cleanedfinaldata/bacterial_refs/Staphylococcus_epidermidis_ASM609437v1_bac.fasta']

print(sw_overfile(reads, hum, contam))
print(round(start_time - time.time(), 3))


