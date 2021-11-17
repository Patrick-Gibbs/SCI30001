# SCI30001 Research Project

Here I include code from my Univerity of Melborune Research Project – Investigating The Potential of Population Specific Homing Drives in Lolium Rigidum.

I pressent two pipelines which pertain to experiment 1 and exerpeiment 2 decribed in the research report.


**Experiment 1**: Calculate the D value for nucleotides in CRISPR/cas9 annealing site for a population.

Dependencies: Python3.9 with Pandas 1.2.4, Biopython  1.79, Numpy 1.19.5.

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



Execution: run Experiment1 with the following: python3 Experiment1.py syncfile.sync 0 

argument 1 is the sync file anylised, and argument two is the population targetted in the sync file (counting from 0 to n-1 in the columns of the sync file).

If a gene annotation file is avialable the follow command will add the gene annotation to the outputted CSVs

