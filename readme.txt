A naive tool for identifying potential DnaA boxes in bacterial genomes.

The tool functions by identifying a genome index of minimum skew. In this case minimum skew is identified by locating the index in the genome at 
which the relative amount of G nucleotides to C nucleotides is minimized. This index is taken to be close to the origin of replication.

Once the ori is located, the tool will search for repeated sequences within a user specified window length. Sequences from starting, ending, and centered around the 
theorized ori index will be searched for the repeated kmers.

The tool accepts certain arguments to make DnaA box searching customizable. It takes an argument to specify the length of the DnaA box, generally 9 bps. It accepts a minimum number of
repeated kmers that indicate a match. It also accepts a hamming distance at which kmers are considered repeats. For instance at a hamming distance of 1, AAT and AAA are considered repeats.

The tool may be utilized via the command line with the following required arguments:
-p, --path = the path to the text file containing the genome in question
-w, --window-length = the length of the window in which to search for repeated sequences
-n, --num-kmers = the minimum number of repeats necessary for a sequence to be considered a potential DnaA box
-d, --dnaa-box-length = the length of DnaA boxes for which to search 
-h, --hamming-distance = the hamming distance at which to consider two sequences to be repeats 

Ex. python find_dnaa_box.py -p Salmonella_enterica.txt -w 1000 -n 6 -d 9 -h 1

The tool prints the theorized DnaA box sequences to the console. 
