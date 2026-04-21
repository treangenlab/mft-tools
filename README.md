# Min-Frame Transformation (MFT)

This code base contains utility functions and analyses for the Min-Frame Transformation (MFT). 

## Overview of Min-Frame Transformation
The MFT allows local transformation of a nucleotide sequence to a character sequence over a separate defined alphabet. This transformation allows effectively masks a large percentage of single-nucleotide mutations and can lead to increased sensitivity when using full-text indexing methods (MUMs, BWT, etc). Downstream this can lead to better accuracy and performance on genome alignment tasks, and potentially more. 

The MFT is similar to the use of minimizers in many ways, however there are some key differences:
1. A value is chosen for every single window, even if it is the same as the neighboring window. This allows for a 1:1 transformation in sequence length
2. The ordering of the kmers is not random and must be defined a-priori
3. Kmers map directly to a reduced alphabet (typically 12-64 characters), rather than a single hash per kmer. This allows for mutations to be masked
