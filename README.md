# SCI30001 Research Project

Here I include code from my Univerity of Melborune Research Project – Investigating The Potential of Population Specific Homing Drives in Lolium Rigidum.

I pressent two pipelines which pertain to experiment 1 and exerpeiment 2 decribed in the research report.



**Dependencies: Python3.9 with Pandas 1.2.4, Biopython  1.79, Numpy 1.19.5.**

**Experiment 1**: Calculate the D value for nucleotides in CRISPR/cas9 annealing site for a population.


Input: The program takes a syncronised (sync) file.

Output: the program outputs a csv called "Results" which contains the following infomation for each nucleotide possition in all CRISPR/cas9 annealing site of the target population:
* scafoldID – The Chromsome or scaffold of the nucelotide possition
* possition – The possition of the CRISPR/cas9 annealing site in the scafold (as specified in the sync file) 
* exact_possition – The possition of the nucleotide in the scafold (as specified in the sync file)
* polymorphism_index – The possition of the nucleotide in the CRISPR/cas9 annealing site (counting from 0)
* target_freq – The frequency of the nucleotide in the target population
* in_gene – spefifies if the nucleotide lays in a (provided anotation file is provided(
* spec_score – The D score for that nucleotide possition
* target_population: restates the population targeted
* 1,2,3...n – gives the frequencey of the most common nucleotide in the target population, in all other populations.

Execution: to run Experiment1 navigate into /Experiment1/ then use the following: python3 Experiment1.py [syncfile.sync] [population targeted] 

argument 1 is the sync file anylised, and argument two is the population targetted in the sync file (counting from 0 to n-1 in the columns of the sync file).

Example Execultion with sample data: nagivagate to /Experiment1/ then run: 

```python3 findDvalue.py ../Sample_Data/Lolium_small.sync 0```

this will give an output Results_0.csv

If a gene annotation file is avialable the running gene_ont.py will add the gene ontolgy row to the output csv. To run it, edit the code such the "file" pertains to a list of integers corrosponding to the populations that need to be converted. for example if one would like to convert population 0 from the example above, ensure file is set to [0]. navigate to /Experiment1/ and use the following command:

```python3 gene_ont.py [annotation file]```

**Experiment 2:**
Code for experment two compares two populations, and finds polymorphic CRISPR/Cas9 annealing sites between them. In the final report for the research project the number of polymorphic CRISPR/Cas9 annealing sites was used. It should be noted that code used here has been readpated serveral times, and thus does not provide the most elequant implementation.

Input: A sync file and two populations (given as index in sync file [0-...])
Output: Two fasta files, one named 30mers.fa and another 23mers.fa as well a file Frequenceys.csv.

23mers.fa gives the sequences of all polymorhic CRISPR annelaing sites (20nt gRNA annealing site + NGG PAM) between the two populations, in the fasta header is the populations number and, whether the sequence is found on the plus or minus strand, and the index in the file to which the annealing site pertains.

30mers.fa gives the same sequence as 23mer plus additional surrounding nucleotides. This file is created such that it can be read by CRISPRon, a program that estimates cleavage effeiency for a CRISPR/cas9 annealing site based on sequence and can be found here https://github.com/RTH-tools/crispron. As it was not mentioned in the research report implmentation of crispron is not included here, but can be provided apon request.

frequenceys.csv gives additional infomation about the annealing sites in 23mers.fa/30mers.fa:
* possition - the genomic location in the given scafold of the 
* polymorphims – the possitions in each CRISPR/cas9 annealing contain polymorphics
* polymorphim_frequency – the frequencey of the nucleotide possitions described in polymorphisms

Execulation – navigate to /Expermiment2/ and use the following command:

```python3 polymorphicSites.py [target population] [off target populations] [sync file]```

target and off target populations are given by integers pertaining to the index in the sync file. For example if one wanted to compared population 0 to population 1 they would use:

```python3 polymorphicSites.py 0 1 ../Sample_Data/Lolium_small.sync```

for doing pair wise comprasons as was conducted in the research report it is recommneded to use GNU parallel. In particualr the following command was used in the research report.

```parallel python3 ./polymorphicSites.py {= 'if($arg[1]==$arg[2]) { skip() }' =} $sync_file ::: {0..63} ::: {0..63}```

by defult results are stored in the 'outputs' directory

**Other files:**

Here I include the distance matrix and fst matrix which were computed to compair pairs of populatoins. each row number pertains to the index provided in the sync file. the first row for both of these files gives the index and was used such that the matrix could be simply imported to pandas, however it can otherwise be ignored.



