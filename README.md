# motif-mark

## About
This oriented object programming script searches for motif sites in gene sequences using a FASTA file as input. The sequences in the input file can contain up to 1000 bases. The exons in the sequence must be upper case, while the intron has to be lower case. The script can search for up to 5 motifs (maximum 10 bases each) and their sequences have to be input as a separated .txt file where each line contains a motif. In addition, the script was written to account for motifs with ambiguous nucleotides (https://en.wikipedia.org/wiki/Nucleic_acid_notation).

If the sequences in the input FASTA file are split along multiple lines, the script will output a file with one-line sequences. Along with that, the output of the script will be a .png image showing the position of the exon and the motifs along the gene with the sizes to scale. Each motif will have a different color, described in the figure legend, and the gene will be identified with a header containing the gene name, chromosome, and location (start and end positions).


## How to use
In order to work, the following dependencies must be installed:

```
Python 3.10.9
Pycairo 1.21.0
```

To run the `motif-mark-oop.py` script you have to include the input FASTA file and the motif file:

``` 
motif-mark-oop.py -f <input_file.fasta> -m <motifs.txt>
```






