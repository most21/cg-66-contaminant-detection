import os
import sys

sys.path.insert(1, '/Users/eleanorhilgart/github/cg-66-contaminant-detection/data_utils')

from data_utils import read_fasta_files
from fasta import FASTA
from fasta import _read_file

def write_solution(ID, seq, out_fh, per_line=60): #took the FASTA writing function from HW4
  offset = 0
  out_fh.write(f">{ID}\n")
  while offset < len(seq):
    line = seq[offset:offset + per_line]
    offset += per_line
    out_fh.write(line + "\n")

'''
This will only need to be used once. The human genome files have a lot of Ns in them, which result from
gaps/inconsistency in the assembly. In order to create in silico fastq reads using the human genomes,
I'm removing the Ns to just get the sequences.
'''

ls = [] #list of input files to be modified
outlst = [] #list of output files to be written to

directory = '/Users/eleanorhilgart/Desktop/CHM1.1_V1_chromosomes_hum_iterate'
for filename in os.listdir(directory):
    if filename.startswith("chr"):
        ls.append(os.path.join(directory, filename))
    elif filename.startswith("noN"):
        outlst.append(os.path.join(directory, filename))
    else:
        continue
ls = sorted(ls)
outlst = sorted(outlst)

for i in range(len(ls)):
  fastafile = read_fasta_files(ls[i]) #read each file
  tup = _read_file(fastafile) #parse the file
  for char in 'Nn':
      tup[0] = tup[0].replace(char, '')
  write_solution(tup[1], tup[0], outlst[i]) #write to the output file
