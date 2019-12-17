# misc-sequence-analysis
Sharing some functions I've written to help me work with genetic sequence data. Most need the "ape" and the "adephylo" package.

## R_concat_func.R
A function that reads in a list of alignments (filenames, alignments must be in PHYLIP), and concatenates the listed files, in the order given, to produce one large concatenated gene alignment. *It requires loading the ape R package.* The version of the function that reads in FASTA as opposed to PHYLIP sequence data is broken, but I've included it in case someone might want to look at it. I tried to be clear as to what the problem is (in the comments).
