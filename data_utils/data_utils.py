from fasta import FASTA
from fastq import FASTQ

def read_fasta_files(files):
    """
        Given a list of FASTA files, return a list of FASTA dataset objects (one object per file).
        Note: this function can also accept a single file as a string, in which case a single FASTA object will be returned.
    """
    if isinstance(files, list):
        return [FASTA(f) for f in files]
    elif isinstance(files, str):
        return FASTA(files)
    raise ValueError(f"Invalid argument type of {type(files)}. Must be str or list")

def read_fastq_files(files):
    """
        Given a list of FASTQ files, return a list of FASTQ dataset objects (one object per file).
        Note: this function can also accept a single file as a string, in which case a single FASTA object will be returned.
    """
    if isinstance(files, list):
        return [FASTQ(f) for f in files]
    elif isinstance(files, str):
        return FASTQ(files)
    raise ValueError(f"Invalid argument type of {type(files)}. Must be str or list")



if __name__ == "__main__":
    ## FASTA Tutorial ##

    # We can create a FASTA object for a particular file
    # Note: We also can pass in a list of FASTA file names and get back a list of FASTA objects
    f = "test_data/fasta_test2.fa" # use whatever FASTA file you like to test this out
    FA = read_fasta_files(f)
    print(FA) # printing a FASTA object displays a nice little summary of the data

    # We can get the number of fragments in the FASTA file using the built-in len()
    length = len(FA)
    print("\nLength = ", length)

    # FASTA objects support efficient, 0-based indexing. FA[frag_num][idx] returns a character
    frag1 = FA[0] # returns a string containing the first fragment in the data
    first_base = FA[0][0] # returns the first character in the first fragment of the FASTA file
    print("First fragment = ", frag1)
    print("First base = ", first_base)

    # We can loop through the FASTA using len() and indexing
    for i in range(len(FA)): # loop over fragments
        for j in range(len(FA[i])): # loop over bases
            base = FA[i][j]
            print(i, j, base)
            break
        break
    print()

    # FASTA objects are also iterable, so you can easily loop through them without knowing their internal representation
    for i, frag in enumerate(FA):
        for j, base in enumerate(frag):
            print(i, j, base)
            break
        break
    print()

    # You can get the list of fragments via the get method
    genome = FA.get_data()
    print(genome)

    print("\n----------------------\n")

    ## FASTQ Tutorial ##

    # FASTQ objects are very similar to FASTA objects, but are slightly more complex.
    # As before, we can create a FASTQ object with a single file or a list of files
    f = "test_data/phix_reads20.fastq" # use whatever FASTQ file you like to test this out
    FQ = read_fastq_files(f)
    print(FQ) # printing displays a summary of the dataset

    # FASTQ objects support len(), which returns the number of reads in the file
    num_reads = len(FQ)
    print("\nNumber of reads = ", num_reads)

    # FASTQ objects support integer indexing.
    # Accesses return a dictionary containing the read, it's id, and its quality values
    first_read_by_idx = FQ[0]
    print(first_read_by_idx)

    # FASTQ objects can also be indexed by read id.
    # Note that this is an O(n) operation, so it is not recommended. If necessary, I'll improve it.
    first_read_by_read_id = FQ["SRR353630.8901"] # this will fail if you use a file other than the phix_reads20.fastq
    print(first_read_by_read_id)

    # FASTQ objects support iteration. The default iterator returns dictionaries as if indexed by integers
    for i, read in enumerate(FQ):
        print(i, read)
        break

    # For convenience, you can also iterate over specific series (read ids, seqs, qualities)
    # To do so, you specify the series as either "id", "seq", or "quality"
    for i, q in enumerate(FQ("id")):
        print(i, q)
        break

    # If you want to access specific series, you can do that too via the get methods
    read_ids = FQ.get_read_sequences()
    print(read_ids[0:5])
