import argparse
import os
# import sys

from data_utils.data_utils import read_fasta_files, read_fastq_files

from kmer import kmer_engine
from fm import fm_engine
from smithwaterman import sw_engine
from minhash import minhash_engine


def parse_args():
    """ Parse command line arguments. """
    parser = argparse.ArgumentParser(description="This is the main driver file for our contaminant detection pipeline.")

    parser.add_argument("--cont-ref", nargs="+", dest="cont_ref", type=str, help="Paths to FASTA files to use as contaminant reference genomes, or a path to a folder containing the FASTA files. Note that folders are not recursively explored.")
    parser.add_argument("--des-ref", nargs="+", dest="des_ref", type=str, help="Paths to FASTA files to use as desired reference genomes (e.g. human), or a path to a folder containing the FASTA files. Note that folders are not recursively explored.")
    parser.add_argument("--query", nargs=1, dest="query", type=str, required=True, help="Path of FASTQ files to check for contamination, or a path to a folder containing the FASTQ files. Only 1 file can be provided at a time.")
    parser.add_argument("--engine", type=str, required=True, choices=["kmer", "fm", "sw", "minhash"], help="Method to use for contaminant detection. Options: kmer = K-mer index, fm = FM index, sw = Smith-Waterman, minhash = Minhash")
    parser.add_argument("--save", dest="save", action="store_true", help="Save cleaned FASTQ file to disk. If omitted, file will not be saved.")

    parser.set_defaults(save=False)

    return parser.parse_args()

def get_fasta_files_from_arg(args):
    """ Parse command line input for ref/query to handle file/folder inputs. """
    files = []
    for arg in args: # for each arg listed (some may be files, others may be folders)
        if os.path.isdir(arg): # if folder, explore it
            for file in os.listdir(arg): # for each item in the folder
                if os.path.isfile(arg + "/" + file) and file != ".DS_Store": # only add it if it's a file (no sub-dirs)
                    files.append(arg + "/" + file)
        elif os.path.isfile(arg) and arg != ".DS_Store":
            files.append(arg)
        else:
            raise IOError("Argument is not a file or directory.")
    return files

# def get_size(obj, seen=None):
#     """
#         Recursively finds size of objects.
#         This function is copied verbatim from the starter code in HW3 in EN.601.446/646.Sp21
#     """
#     size = sys.getsizeof(obj)
#     if seen is None:
#         seen = set()
#     obj_id = id(obj)
#     if obj_id in seen:
#         return 0
#     # Important mark as seen *before* entering recursion to gracefully handle
#     # self-referential objects
#     seen.add(obj_id)
#     if isinstance(obj, dict):
#         size += sum([get_size(v, seen) for v in obj.values()])
#         size += sum([get_size(k, seen) for k in obj.keys()])
#     elif hasattr(obj, '__dict__'):
#         size += get_size(obj.__dict__, seen)
#     elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
#         size += sum([get_size(i, seen) for i in obj])
#     return size

def write_clean_fastq_file(filename, results):
    """ Write a FASTQ file without the contamintated reads. """
    filename_no_ext = "".join(filename.split(".")[:-1])
    with open(filename_no_ext + "_clean.fastq", "w") as f:
        output = ""
        # Iterate through non-contaminated/unassigned reads
        for read in results["Desired"]:
            r_id, r_seq, r_qual = read["id"], read["seq"], read["quality"]
            output += f"@{r_id}\n{r_seq}\n+\n{r_qual}\n"
        for read in results["Unassigned"]:
            r_id, r_seq, r_qual = read["id"], read["seq"], read["quality"]
            output += f"@{r_id}\n{r_seq}\n+\n{r_qual}\n"

        # Write to file
        f.write(output)


def main():
    # Parse command line args
    args = parse_args()

    # Get all files specified in the command line
    cont_fasta_files = get_fasta_files_from_arg(args.cont_ref)
    des_fasta_files = get_fasta_files_from_arg(args.des_ref)
    fastq_file = args.query[0] # even though it's just 1 file, it comes in as a list of size 1

    # Read all reference genome FASTA files and all query FASTQ files
    print("Creating FASTA objects for all reference files...")
    cont_FA = read_fasta_files(cont_fasta_files)
    des_FA = read_fasta_files(des_fasta_files)
    print("Creating FASTQ object for query file...")
    FQ = read_fastq_files(fastq_file)

    # Run contaminant detection method based on the selected engine
    engine = args.engine
    if engine == "kmer":
        print("Running k-mer index engine...\n")
        results = kmer_engine(FQ, des_FA, cont_FA)
    elif engine == "fm":
        results = fm_engine(FQ) # TODO: args
        #raise NotImplementedError("TODO: Implement FM index engine")
    elif engine == "sw":
        print("Running Smith-Waterman engine...\n")
        results = sw_engine(FQ, des_FA, cont_FA)
    elif engine == "minhash":
        print("Running MinHash engine...\n")
        results = minhash_engine(cont_FA, des_FA, FQ)
    else:
        # The command line args should already validate the engine selection, so this line will likely never be run.
        raise ValueError("Invalid engine for contaminant detection.")

    # Display results
    contam_percent = round(100 * len(results["Contaminated"]) / len(FQ), 2)
    des_percent = round(100 * len(results["Desired"]) / len(FQ), 2)
    unassigned_percent = round(100 * len(results["Unassigned"]) / len(FQ), 2)

    print(f"Contaminated: {contam_percent}%")
    print(f"Desired: {des_percent}%")
    print(f"Unassigned: {unassigned_percent}%")

    # Optionally create output file with contaminated reads removed
    if args.save:
        print("\nCreating cleaned file...")
        write_clean_fastq_file(fastq_file, results)


if __name__ == "__main__": main()
