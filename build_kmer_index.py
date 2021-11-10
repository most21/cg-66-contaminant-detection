#Build a basic k-mer index from a fasta file

#import desired functions from data_utils
from data_utils import read_fasta_files

def build_kmer_index(fasta_name, k):
    """ Construct a basic kmer index from an input FASTA file and a kmer length. """
    #import the desired fasta file and make a fasta object using the data utils function
    FA = read_fasta_files(fasta_name)

    #make an empty dictionary for the build_kmer_index
    #keys will be kmers and entries will be a list of where in the fasta the kmer starts
    idx = {}
    for i in range(len(FA)-k+1):
        mer = FA[i:i+k]
        #make new entries for new kmers
        if mer not in idx:
            idx[mer] = [i]
        #append the location of repeated kmers to the list
        elif mer in idx:
            idx[mer].append(i)

    #return the completed kmer index
    return idx
