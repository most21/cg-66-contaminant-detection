# cg-66-contaminant-detection
Computational Genomics Final Project - Group 66 - Contaminant Detection  
Team Members: Eleanor Hilgart, Matthew Ost, Gabe Shatkin, & Alexey Solganik

## Install + Run
- Install requirements with pip. The only requirements are Python 3.7+, numpy, xxhash, and [filprofiler](https://pythonspeed.com/fil/docs/index.html).
```shell
pip install -r requirements.txt
```

- Run `main.py` with the proper command line args. Do `python main.py -h` for more information.
```shell
python main.py ...
```

- To run on a very small example, first unzip the tinydataexample.zip folder. Then specify an engine (kmer, fm, sw, minhash) and run:
```shell
unzip tinydataexample.zip
python main.py --cont-ref tinydataexample/Mfermentansbac1_cut.fasta --des-ref tinydataexample/noN_chr1_cut.fasta --query tinydataexample/chr1_Mfermentans_8020_hiseq_reads_tinycut_R1.fastq --engine <TODO: ENGINE>
```

## Data
Since we have about 12 GB of data, we are storing all the files in Google Drive at this [link](https://drive.google.com/drive/folders/1MQlJKj1cV6ziyyoWNRC5i00vMFAE5WlQ?usp=sharing).  
A document explaining the data can be found in the drive and is also linked [here](https://tinyurl.com/3kmn8k3r).  
This repository includes a zip folder called tinydataexample.zip, which contains a very small set of files for testing.
All our methods can run in a reasonable amount of time on this data.

## Profiling
Our paper displays a table containing time and memory measurements for each engine on the tiny data example shown above. To get the memory results for each engine, run the line below with the engine of your choice:

```shell
python -m filprofiler run main.py --cont-ref tinydataexample/Mfermentansbac1_cut.fasta --des-ref tinydataexample/noN_chr1_cut.fasta --query tinydataexample/chr1_Mfermentans_8020_hiseq_reads_tinycut_R1.fastq --engine <TODO: engine>
```

To get the time results for each engine, we used the Unix time tool and reported the "real" time output.
```shell
time python main.py --cont-ref tinydataexample/Mfermentansbac1_cut.fasta --des-ref tinydataexample/noN_chr1_cut.fasta --query tinydataexample/chr1_Mfermentans_8020_hiseq_reads_tinycut_R1.fastq --engine <TODO: engine>
```

## Repository Structure
- `main.py` is the driver script and entry point to our software
- `kmer.py` contains the code for the k-mer index method
- `fm.py` contains the code for the FM index method
- `smithwaterman.py` contains the code for the Smith-Waterman method
- `minhash.py` contains the code for the MinHash method
- `data_utils` is a directory containing the supporting tools for working with data
  - `data_utils/data_utils.py` contains a few key functions to create FASTA/Q objects and load/save Python objects from/to disk
  - `data_utils/dataset.py` contains an abstract base class for Dataset objects
  - `data_utils/fasta.py` contains the FASTA class definition
  - `data_utils/fastq.py` contains the FASTA class definition
  - `data_utils/data_removeNs.py` contains the code to remove no-confidence bases from reads
- `cache/minhash` is a directory containing MinHash sketches for many of the reference genomes in our dataset
- `tinydataexample.zip` contains small data files for a working example of our code

A user only needs to download the desired data and run `main.py` with the proper command line arguments to run each of our methods.
