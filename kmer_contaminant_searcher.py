from data_utils import read_fastq_files
from build_kmer_index import build_kmer_index

def kmer_cont_search(fq_name, des_ref, cont_refs, k):
    #at some point want to include input for approximate match threshold
    """Count the reads from a fastq object that allign to the desired index or contaminant indexes"""

    #make fastq object
    FQ = read_fastq_files(fq_name)
    #make index from desired reference 
    des_idx = build_kmer_index(des_ref, k)
    #make indexes from contamination references
    cont_idxs = []
        #check if we have a list of contamination references
    if type(cont_refs) is list:
        for cont_ref_name in cont_refs:
            cont_idxs.append(build_kmer_index(cont_ref_name, k))
            #check if we have a single contamination reference
    if type(cont_refs) is str:
        cont_idxs.append(build_kmer_index(cont_refs, k))

    #set the counts of desired and contaminant reads
    des_count = 0
    cont_count = 0
    unas_count = 0

    #start with exact match, but want to consider approximate later

    #look at each read in the fastq file. 
    #Using kmer method, try to align it to the desired and contaminant references
    for read in FQ:
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
        #if all kmers are in the reference, check if they are all next to each other
        if len(kmer_locs) == len(read_kmers) and len(kmer_locs) > 0:
            for start in kmer_locs[0]:
                is_in = 1
                check = start + 1
                for i in range(1,len(kmer_locs)):
                    if check in kmer_locs[i]:
                        check = check + 1 
                    elif check not in kmer_locs[i]:
                        is_in = 0    
                        break
                #if all looks good, increase the desired reference alligned reads count
                if is_in == 1:
                    des_count = des_count + 1
                    in_des = 1
                    break
                
        #check if the k_mers are in the contamination references
        kmer_locs = []
        for mer in read_kmers:
            for idx in cont_idxs:
                if mer in idx:
                    kmer_locs.append(idx[mer])
                #if all kmers are in the reference, check if they are all next to each other
            ###check back at the indent spacing here
            if len(kmer_locs) == len(read_kmers) and len(kmer_locs) > 0:
                for start in kmer_locs[0]:
                    is_in = 1
                    check = start + 1
                    for i in range(1,len(kmer_locs)):
                        if check in kmer_locs[i]:
                            check = check + 1 
                        elif check not in kmer_locs[i]:
                            is_in = 0       
                            break
                    if is_in == 1:
                        cont_count = cont_count + 1
                        in_cont = 1
                        break
    
        #check if the read was unassigned to either a desired or contamination
        if in_des == 0 and in_cont == 0:
            unas_count = unas_count + 1

    #return a list of the counts of reads mapped to the desired reference, a contamination reference, or unassigned
    results = [des_count, cont_count, unas_count]
    return(results)
            
