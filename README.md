# misc-sequence-analysis
Sharing some functions I've written to help me work with genetic sequence data. Most need the "ape" or the "adephylo" package.

## concat_func.py
A function that uses Biopython to convert a list of alignment filenames into one AlignIO object. I've avoided making it into a standalone script to allow for more flexibility (ex. the version here reads in alignments in the FASTA format, but AlignIO can read a large diversity of alignment formats - see (here)[https://biopython.org/wiki/AlignIO])

## R_concat_func.R
A function that reads in a list of alignments (filenames, alignments must be in PHYLIP), and concatenates the listed files, in the order given, to produce one large concatenated gene alignment. *It requires loading the ape R package.* The version of the function that reads in FASTA as opposed to PHYLIP sequence data is broken, but I've included it in case someone might want to look at it. I tried to be clear as to what the problem is (in the comments). This function is outdated and I would recommend not using it (2025)
