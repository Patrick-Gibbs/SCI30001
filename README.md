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

Execution: to run Experiment1 navigate into /Experiment1/ then use the following: python3 Experiment1.py [syncfile.sync] [population targeted] 

argument 1 is the sync file anylised, and argument two is the population targetted in the sync file (counting from 0 to n-1 in the columns of the sync file).


Example Execultion with sample data: nagivagate to /Experiment1/ then run: python3 findDvalue.py ../Sample_Data/Lolium_small.sync 0
this will give an output Results_0.csv

If a gene annotation file is avialable the running gene_ont.py will add the gene ontolgy row to the output csv. To run it, edit the code such the "file" pertains to a list of integers corrosponding to the populations that need to be converted. for example if one would like to convert population 0 from the example above, ensure file is set to [0]. navigate to /Experiment1/ and use the following command:

python3 gene_ont.py [annotation file]







