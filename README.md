# cg-66-contaminant-detection
Computational Genomics Final Project - Group 66 - Contaminant Detection

## Install + Run
- Install requirements. The only requirements are Python 3.7+, numpy, and xxhash.
```shell
pip install -r requirements.txt
```

- Run `main.py` with the proper command line args. Do `python main.py -h` for more information.
```shell
python main.py ...
```

- To run on a very small example, specify an engine (kmer, fm, sw, minhash) and run:
```shell
python main.py --cont-ref tinydataexample/Mfermentansbac1_cut.fasta --des-ref tinydataexample/noN_chr1_cut.fasta --query tinydataexample/chr1_Mfermentans_8020_hiseq_reads_tinycut_R1.fastq --engine <TODO: ENGINE>
```
