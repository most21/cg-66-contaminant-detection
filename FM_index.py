#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 17:05:03 2021

@author: Alex
"""

import copy as cp
import numpy as np
import math

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



######################################
#Functions I wrote

def totsBWT(bw):
    '''Map of characters to number of times it appears '''
    tots = dict()
    for c in bw:
        if c not in tots:
            tots[c] = 0
        tots[c] += 1
    return tots


def firstColMod(tots):
    ''' Return map from character to the range of rows prefixed by
        the character. I changed the function to 1-indexed vs 
        0-indexed before'''
    first = {}
    totc = 1
    for c, count in sorted(tots.items()):
        first[c] = (totc, totc + count)
        totc += count
    return first


def modifiedRankBWT(bw,step):
    tots = totsBWT(b)
    f = firstColMod(tots)
    tots = dict()
    total = 0

    for i in f:
        f[i] = []
    f1 = cp.deepcopy(f)
    f2 = cp.deepcopy(f)
            
    first = True
    for i in range(0,len(bw)):
        c = bw[i]
        if bw[i] not in tots:
            tots[c] = 0
        
        total += 1
        tots[c] += 1
        f1[c].append([total, tots[c]])
        
        
        for j in f1:
            if j != c:
                #print(first)
                if len(f1[j]) == 0:
                    if first == False:
                        last_entry_from_main_table = cp.deepcopy(f[j][-1])
                        last_entry_from_main_table[0] +=1
                        f1[j].append(last_entry_from_main_table)               
                    else:
                       
                       f1[j].append([1,0])
                else:
                    previous_entry =  cp.deepcopy(f1[j][-1])
                    previous_entry[0] +=1
                    f1[j].append(previous_entry)  
        first = False            
        if i % step == 0:
            for k in f:
                f[k].append(f1[k][-1])
            f1 = cp.deepcopy(f2)  
    for i in f:
        f[i] = np.asarray(f[i])
       
    return f

t = 'abaabab$'
b = bwtViaBwm(t)

tots = totsBWT(b)
#print(ranks)

print(b)

print(firstColMod(tots))



def find_nearest(array,value):
    idx = np.searchsorted(array[:,0], value, side="left")
    
    if idx > 0 and (idx == len(array[:,0]) or math.fabs(value - array[:,0][idx-1]) < math.fabs(value - array[:,0][idx])):
        return array[idx-1]
    else:
        return array[idx]
    
    


def queryBWT(pattern, index_table, bw):
    
    # Get the last char of the pattern
    """ Reverse the order of the string for ease"""
    reverse_pattern = pattern[::-1]
    length = len(reverse_pattern)
    
    start = 0
    """ Start with the 1st char of theis string"""
    starting_char = reverse_pattern[start]
    
    tots = totsBWT(b)
    f = firstColMod(tots)
    
    total = 0
    """ Get the first range of characters from the first column"""
    first_range = f[starting_char]
    #print(first_range)
    final_rank = 0
   
    next_range = range(first_range[0], first_range[1])       
    next_range1 = []
    #print(next_range)
    for k in range(1, len(pattern)):
        start += 1
        """ If there is a next character in the reverse_pattern"""
        if start < length:
            """ Go through the range of characters in the first column"""
            for j in range(len(next_range)): 
                # bw is 0-indexed so have to subtract 1
                """ Goet the next character"""
                c = reverse_pattern[start]
                i = next_range[j]
                
                
                """ Compare the next character with a character at the same index 
                from BWT. Have to subtract 1 because BWT is 0-indexed."""
                if bw[i-1] == c:
                    
                    """In the dictonary go to the 'column' of the desired 
                    character and find the nearest available total index"""
                    closest_index = find_nearest(index_table[c], i)
                    """If the nearest available total index is the actual total
                    index then we can return the corresponding char-specific index.
                    That will be the actual index of the character."""
                    if closest_index[0] == i:
                        final_rank = closest_index[1]
                    """If the nearest available total index is smaller than 
                    the actual total index, we have to walk down bwt"""
                    
                    if closest_index[0] < i:
                        # Subtract 1 because bw is 0-indexed
                        closest_in_bw = closest_index[0] 
                        final_rank = closest_index[1]
                        for j in range(closest_in_bw, i):
                            if bw[j] == c:
                                final_rank += 1       
                    if closest_index[0] > i:
                        # Subtract 1 because bw is 0-indexed
                        closest_in_bw = closest_index[0] 
                        final_rank = closest_index[1]
                        for j in range(i, closest_in_bw):
                            if bw[j] == c:
                                final_rank -= 1    
                
                    """After we have found the index of the character, need to 
                    find where it is in the left column"""
                
                    transition_index = f[c][0] + final_rank -1
                    next_range1.append(transition_index) 
                
            next_range = cp.deepcopy(next_range1)        
            next_range1 = []    
            
            
    return len(next_range)
           
        
            
        
    
    
   
m = modifiedRankBWT(b,3)

print(m)

print(queryBWT("abaa", m, b))





#def queryFM(bw):
    #''' Given a bw and another string check if 


#print(b)
#print(reverseBwt(bwtViaBwm("Tomorrow_and_tomorrow_and_tomorrow$")))

