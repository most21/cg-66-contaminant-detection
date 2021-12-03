from data_utils import read_fastq_files
from build_kmer_index import build_kmer_index

def kmer_cont_search(fq_name, des_ref, cont_refs, k, tolerance):
    #at some point want to include input for approximate match threshold
    """Returns three lists containing which (one-indexed) reads allign to the desired reference, contaminant reference, or neither"""
   

    #make fastq object
    FQ = read_fastq_files(fq_name)
    #make index from desired reference 
    otpt = build_kmer_index(des_ref, k)
    des_idx = otpt[0]
    des_FA = otpt[1]
    #make indexes from contamination references
    cont_idxs = []
    cont_FA = []
        #check if we have a list of contamination references
    if type(cont_refs) is list:
        for cont_ref_name in cont_refs:
            otpt = build_kmer_index(cont_ref_name, k)
            cont_idxs.append(otpt[0])
            cont_FA.append(otpt[1])
            #check if we have a single contamination reference
    if type(cont_refs) is str:
        otpt = build_kmer_index(cont_refs, k)
        cont_idxs.append(otpt[0])
        cont_FA.append(otpt[1])
        
        
    #set the counts of desired and contaminant reads
    des_count = 0
    cont_count = 0
    unas_count = 0
    
    
    #start with exact match, but want to consider approximate later

    #look at each read in the fastq file. 
    #Using kmer method, try to align it to the desired and contaminant references
    which_read = 0 #keep track of which read we're on
    unassigned_reads = [] #keep track of which reads have been assigned to what
    good_reads = []
    cont_reads = []
    for read in FQ:
        which_read += 1
        unassigned_reads.append(which_read)
        #keep track of where each read has been assigned
        in_des = 0
        in_cont = 0
        read_kmers = []
        #make kmers from the read
        for i in range(len(read)-k+1):
            read_kmers.append(read[i:i+k])

        #check if the k_mers are in the desired reference
        kmer_locs = []
        for mer in read_kmers:
            if mer in des_idx:
                kmer_locs.append(des_idx[mer])
                
        #Check for alignments to the desired reference 
        if len(kmer_locs) >= len(read_kmers)-((k-1)*tolerance) and len(kmer_locs) > 0:
            #need to check the reference and calculate the hamming distance each start, as well as k, 2k, ... and (tolerance)k before it
            for start in kmer_locs[0]:
                starts = [start]
                for i in range(tolerance):
                    if starts[len(starts)-1] - k >= 0:
                        starts.append(starts[len(starts)-1] - k)
                        
                #need to compare the read to each read length section of the reference begining at each start point
                for strt in starts:
                    num_mismatch = 0
                    compare = des_FA[strt:strt+len(read)]
                    for base_num in range(len(read)):
                        base = read[base_num]
                        if base != compare[base_num]:
                            num_mismatch += 1
                        if num_mismatch > tolerance:
                            break
                    if num_mismatch <= tolerance and which_read not in good_reads:
                        des_count = des_count + 1
                        in_des = 1
                        good_reads.append(which_read)
                        if which_read in unassigned_reads:
                            unassigned_reads.remove(which_read)
                        break
                
        #check if the k_mers are in the contamination references only if the read isn't in the desired ref
        if in_des != 1:
            kmer_locs = []
            for mer in read_kmers:
                if mer in cont_idxs[0].keys():
                    kmer_locs.append(cont_idxs[0][mer])
            
                #for idx in cont_idxs:
                    #if mer in idx:
                        #kmer_locs.append(idx[mer])
                    #if all kmers are in the reference, check if they are all next to each other
                ###Deleted code here(find in bottom cell)
###need to add in aproximate matching with contminant reference
            if len(kmer_locs) >= len(read_kmers)-((k-1)*tolerance) and len(kmer_locs) > 0: #filter which reads could potentially be a match by pigeonhole
                new_starts = kmer_locs[0]
                    #need to check the reference and calculate the hamming distance each start, as well as k, 2k, ... and (tolerance)k before it
                    #need to compare the read to each read length section of the reference begining at each start point
                for strt in new_starts:
                    num_mismatch = 0
                    compare = cont_FA[0][strt:strt+len(read)]
                        
                    for base_num in range(len(read)):
                        base = read[base_num]
                        if base != compare[base_num]:
                            num_mismatch += 1
                        if num_mismatch > tolerance:
                            break
                    if num_mismatch <= tolerance and which_read not in cont_reads:
                        cont_count = cont_count + 1
                        in_cont = 1
                        cont_reads.append(which_read)
                        if which_read in unassigned_reads:
                            unassigned_reads.remove(which_read)
                        break
    
        #check if the read was unassigned to either a desired or contamination
        if in_des == 0 and in_cont == 0:
            unas_count = unas_count + 1

    #return a list of the counts of reads mapped to the desired reference, a contamination reference, or unassigned
    results = [good_reads, cont_reads, unassigned_reads]
    return(results)
            
