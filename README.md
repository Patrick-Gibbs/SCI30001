# SCI30001 Research Project

Here I include code from my Univerity of Melborune Research Project – Investigating The Potential of Population Specific Homing Drives in Lolium Rigidum.

I pressent two pipeline which pertain to experiment 1 and exerpeiment 2 decribed in the research report.


**Experiment 1**: Calculate the D value for nucleotides in CRISPR/cas9 annealing site for a population.

Dependencies: Python3.9 with Pandas 1.2.4, Biopython  1.79, Numpy 1.19.5.

Input: The program takes a syncronised (sync) file.

Output: the program outputs a csv called "Results" which gives the following:
* abc

Execution: run Experiment1 with the following: python3 Experiment1.py syncfile.sync 0 

argument 1 is the sync file anylised, and argument two is the population targetted in the sync file (counting from 0 to n-1 in the columns of the sync file).
