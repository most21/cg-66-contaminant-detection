#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 17:05:03 2021

@author: Alex
"""

import copy as cp

def rotations(t):
    ''' Return list of rotations of input string t '''
    tt = t * 2
    return [ tt[i:i+len(t)] for i in range(0, len(t)) ]

def bwm(t):
    ''' Return lexicographically sorted list of t’s rotations '''
    return sorted(rotations(t))

def bwtViaBwm(t):
    ''' Given T, returns BWT(T) by way of the BWM '''
    return ''.join(map(lambda x: x[-1], bwm(t)))


def rankBwt(bw):
    ''' Given BWT string bw, return parallel list of B-ranks.  Also
        returns tots: map from character to # times it appears. '''
    tots = dict()
    ranks = []
    for c in bw:
        if c not in tots:
            tots[c] = 0
        ranks.append(tots[c])
        tots[c] += 1
    return ranks, tots
   


def firstCol(tots):
    ''' Return map from character to the range of rows prefixed by
        the character. '''
    first = {}
    totc = 0
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count)
        totc += count
    return first


t = 'abaabab$'
b = bwtViaBwm(t)

ranks, tots = rankBwt(b)
print(ranks)

print(b)
print(firstCol(tots))


def reverseBwt(bw):
    ''' Make T from BWT(T) '''
    ranks, tots = rankBwt(bw)
    first = firstCol(tots)
    rowi = 0 # start in first row
    t = '$' # start with rightmost character
    print(bw[rowi])
    while bw[rowi] != '$':
        #print(bw[rowi])
        c = bw[rowi]
        
        t = c + t # prepend to answer
        # jump to row that starts with c of same rank
        rowi = first[c][0] + ranks[rowi]
    return t


def modifiedRankBWT(bw,step):
    _, tots = rankBwt(b)
    f = firstCol(tots)
    tots = dict()

    for i in f:
        f[i] = []
    f1 = cp.deepcopy(f)
    f2 = cp.deepcopy(f)
            
    first = True
    for i in range(0,len(bw)):
        c = bw[i]
        if bw[i] not in tots:
            tots[c] = 0
        f1[c].append(tots[c]+1)
        tots[c] += 1

        for j in f1:
            if j != c:
                #print(first)
                if len(f1[j]) == 0:
                    if first == False:
                        f1[j].append(f[j][-1])               
                    else:
                       f1[j].append(0)
                else:
                    f1[j].append(f[j][-1])  
        first = False            
        if i % step == 0:
            for k in f:
                f[k].append(f1[k][-1])
            f1 = cp.deepcopy(f2)  
       
    return f







#def queryFM(bw):
    #''' Given a bw and another string check if 


#print(b)
#print(reverseBwt(bwtViaBwm("Tomorrow_and_tomorrow_and_tomorrow$")))

