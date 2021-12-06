#from data_utils import read_fastq_files
#from build_kmer_index import build_kmer_index

def build_kmer_index(FA, k):
    """ Construct a basic kmer index from an input FASTA file and a kmer length. """
    #make an empty dictionary for the build_kmer_index
    #keys will be kmers and entries will be a list of where in the fasta the kmer starts
    idx_list = [] #list to store all the indices created from a multi fastsa file
    FA_list = [] #list to seperately store the FA's from a multi fasta file
    for frag in FA: # TODO: I (Matthew) added this to support multi-fasta files (data is a list of strings instead of 1 string)
        idx = {}
        for i in range(len(frag)-k+1):
            mer = frag[i:i+k]
            #make new entries for new kmers
            if mer not in idx:
                idx[mer] = [i]
            #append the location of repeated kmers to the list
            elif mer in idx:
                idx[mer].append(i)
        idx_list.append(idx)
        FA_list.append(frag)

    #return the completed kmer index
    return idx_list, FA_list

def kmer_engine(FQ, des_ref, cont_refs, k = 10, tolerance = 5): #requires the desired reference to be input as a single fasta file, contaminant can be multi fasta
    #Inputs must be fastq and fasta objects, with integers for k and tolerance
    """Returns three dictionaries containing which reads allign to the desired reference, a contaminant reference, or neither"""

    #make indexes from desired references
    otpt = build_kmer_index(des_ref, k)
    des_idxs = otpt[0]
    des_FAs = otpt[1]

    #make indexes from contamination references
    otpt2 = build_kmer_index(cont_refs, k)
    cont_idxs = otpt2[0]
    cont_FA = otpt2[1]

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
    count_good = []
    count_cont = []
    
    for entry in FQ:
        read = entry['seq']
        which_read += 1
        
        #keep track of where each read has been assigned
        in_des = 0
        in_cont = 0
        read_kmers = []
        #make kmers from the read
        for i in range(len(read)-k+1):
            read_kmers.append(read[i:i+k])

        idx_count = 0
        for des_idx in des_idxs:
            #check if the k_mers are in the desired reference
            kmer_locs = []
            for mer in read_kmers:
                if mer in des_idx:
                    kmer_locs.append(des_idx[mer])

            #Check for alignments to the desired reference
            if len(kmer_locs) >= len(read_kmers)-((k)*tolerance) and len(kmer_locs) > 0:
                #need to check the reference and calculate the hamming distance each start, as well as k, 2k, ... and (tolerance)k before it
                for start in kmer_locs[0]:
                    starts = [start]
                    for i in range(tolerance):
                        if starts[len(starts)-1] - k >= 0:
                            starts.append(starts[len(starts)-1] - k)

                    #need to compare the read to each read length section of the reference begining at each start point
                    for strt in starts:
                        num_mismatch = 0
                        compare = des_FA[idx_count][strt:strt+len(read)]
                        for base_num in range(len(read)):
                            base = read[base_num]
                            if base_num >= len(compare):
                                #num_mismatch = tolerance +1
                                break
                            else:
                                if base != compare[base_num]:
                                    num_mismatch += 1
                                if num_mismatch > tolerance:
                                    break
                     #add aligned reads to the dictionary to output
                    if num_mismatch <= tolerance and which_read-1 not in count_good:
                        des_count = des_count + 1
                        in_des = 1
                        count_good.append(which_read-1)
                        good_reads.append(FQ[which_read-1])
                        #if FQ[which_read-1] in unassigned_reads:
                            #unassigned_reads.remove(FQ[which_read-1])
                        break
            idx_count += 1
#check for each contaminant reference genome provided
        count = 0
        for cont_idx in cont_idxs:
            if in_des != 1:
                kmer_locs = []
                for mer in read_kmers:
                    if mer in cont_idx.keys():
                        kmer_locs.append(cont_idx[mer])

            #approximate matching with contminant reference
                if len(kmer_locs) >= len(read_kmers)-((k)*tolerance) and len(kmer_locs) > 0: #filter which reads could potentially be a match by pigeonhole
                    new_starts = kmer_locs[0]
                        #need to check the reference and calculate the hamming distance each start, as well as k, 2k, ... and (tolerance)k before it
                        #need to compare the read to each read length section of the reference begining at each start point
                    for strt in new_starts:
                        num_mismatch = 0
                        compare = cont_FA[count][strt:strt+len(read)]
                        for base_num in range(len(read)):
                            base = read[base_num]
                        #to prevent checking  past edge of reference
                            if base_num >= len(compare):
                                #num_mismatch = tolerance +1
                                break
                            else:
                                if base != compare[base_num]:
                                    num_mismatch += 1
                                if num_mismatch > tolerance:
                                    break

                    #record reads assigned to a contaminant
                    if num_mismatch <= tolerance and which_read-1 not in count_cont:
                        cont_count = cont_count + 1
                        in_cont = 1
                        count_cont.append(which_read-1)
                        cont_reads.append(FQ[which_read-1])
                        #if FQ[which_read-1] in unassigned_reads:
                            #unassigned_reads.remove(FQ[which_read-1])
                        break
            count += 1
        #check if the read was unassigned to either a desired or contamination
        if in_des == 0 and in_cont == 0:
            unassigned_reads.append(FQ[which_read - 1])
            unas_count = unas_count + 1

    #return a list of the counts of reads mapped to the desired reference, a contamination reference, or unassigned
    results = {"Desired": good_reads, "Contaminated": cont_reads, "Unassigned": unassigned_reads}
    return results
