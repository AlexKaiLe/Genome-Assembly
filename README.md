# README

## How to run the program:

There are 2 shell script files titled debruijn.sh and assembly.sh.

### debruijn.sh 
To run the debruijn script, run the script in the terminal `$ sh debruijn.sh <reads.txt>` where reads.txt is a text file containing a list of DNA reads of constant length (separated by newlines). Example: `sh debruijn.sh test_cases/debruijn/reads2.txt`. The shell script calls a python file that reads in k-mers, convert them into the edges of an internal representation of a de Bruijn graph, simplify this internal representation, and finally print out all possible most likely (i.e., maximum-length) target sequences which can be inferred from the input reads. The program also generates a DOT file called debruijn.dot representing the fully simplified de Bruijn graph. The program simplifies the debruijn graph where every “singleton” is merged into one node and the edges in the resulting graph are updated to reflect the full overlap of the final node labels. 

### assembly.sh
To run the assembly script, run the script in the terminal `$ sh assembly.sh <reads.txt> <vector.txt> <contam-k> <k-mer length> <k> <t> <d>`. Example: `sh assembly.sh test_cases/assembly/sampleReads.txt test_cases/assembly/vector.txt 10 30 17 2 2`. The shell script calls a python file that reads in k-mers from reads and decontaminates the reads with a given vector. The program finds the longest contig of these sequences and returns that sequence.

## Workflow:
I took 3 files and combined them in sequential order to find the final inferred sequence of a sequence of various k-mers. I first took the reads files and removed any of the reads which contained contamination. I did this by running all the reads through the contamination.py file, which returns a list of contaminated reads. Once I removed all the contaminated reads, I determined the reads which I needed to correct. Because I do not know the actual length of the fully sequenced sampleRead.txt, I have to take in all the kmers derived from the reads to sequence to find the longest contig. I place decontaminated reads into the correction.py file, which replaces all infrequent kmers with frequent kmers throughout the reads. Now I have a list of decontaminated and corrected reads. Afterward, I generate kmers from this list of a specified length. These kmers are then placed into the debruijn.py file that then generates the longest contiguous segment that can be inferred from these kmers.

I chose to do contamination before correction because I do not want to be correcting a contamination region. If I corrected a contamination region by accident, then our contamination algorithm would not be able to determine if the read was contaminated or not because that region would have changed as a result of the correction. It is possible that a contamination region may have been altered and thus, preventing the algorithm from detecting the contamination. In this case, using correction to rectify this alteration would be helpful. However, this case is more unlikely to occur and it would prove to be more harmful if I accidentally introduced more alterations in the contamination regions that would prevent the algorithm from identifying it and removing it. 

At the end of the correction phase, I pass all the corrected reads into the de Bruijn file to create the graph and find the longest contig that can be inferred. 

For the sampleReads.txt file, I used the parameters [10, 30, 17, 2, 2]. 
- 10 = kmer length to determine seeding in contamination 
- 30 = kmer length of reads being used to correct and construct a De Bruijn graph 
- 17 = kmer length to determine seeding in correction 
- 2 = threshold value to determine the frequency of common kmer values in correction
- 2 = difference value for how much a replaceable kmer can differ by in correction


## Detecting Errors by optimizing variables for k (k-mer length) and t (threshold for frequent k-mers)

Analyze the hiv-1_genome.txt that contains 9,181 bp RNA genome of HIV-1using unitary.m, which is a DNA scoring matrix. Simulate extraction of 50-mer reads from the HIV genome  where I mutated these reads randomly based on a 1% error rate

P = ​​probability of number of reads covering a given position of the sequence
P=(e^(-a)a^k)/k! where a = N/(LG)
To find the min number of reads covering the entire sequence 
- K = 0 → no reads cover the entire sequence 
- L = 50 → length of every read
- G = 9181 → length of the entire sequence
- N = ? → The number of reads covering the entire length of the sequence 
- 1-P = 1 - the probability that no reads cover the entire sequence = 0.99

Ensure that at least 99.9% of the target sequence is expected to be covered
1 - (e^(-a)a^k)/k!= 0.999
 
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">

- 1 - (e^(-a)a^0)/0! = 0.999
- 1 -e^(-a) = 0.999
- 1 e^(-N/(LG)) = 0.999
- 0.001 =e^(-N/(LG))
- ln(0.001) = -N/(LG)
- (-G • ln(0.001)) / L = N
- (-9181 • ln(0.001)) / 50 = 1268.4 ≈ 1269

Thus, 1269 reads of 50-mers to ensure that 99.9% of the target sequence is expected to be covered

Values of 
- Nmin = 1269
- Sm = 49.981087470449175
- max(Sk) = 49.99211977935382 where k = 22, t = 4
- max(St) = 50 where k = [17-25], t = 2