from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from Population_Sequences import Population_Sequences
from Constants import Constants
from Read import Read
from Thirtymer import Thirtymer
from time import time
import sys

"""
This program finds fixed differences in then CRISPR annealing sites between
a target and one off target population. It then produces the results in a
fasta format which can be read by CRISPRon program to calculate the cleavage
effiencey of these sites.

additaionly a 'frequenceys' file is produced with suplementry data, this can
be combined with the output file from CRISPRon create a csv with more comprehenive
infomation regarding the polymorphic annealing sites found.

"""

#keeps track of how long the program takes to run.
t0 = time()
"""                          -----Step 1-----
    Read in sync file to a refference of the following type:
    {scafold_1: [fixed sequence of population_1, fixed sequence of population_2 ...],
    scafold_2: [fixed sequence of population_1, fixed sequence of population_2 ...],
    ...}
    fixed sequences is the consensus sequence for each population, if a locus
    if not fixed within a popualtion for a given nucleotide (<95%) then '-' is
    put into the sequence.

    sequences are in the biopython seq format
"""
filename = sys.argv[3] # sync file to be read in
other_population = sys.argv[2] # offtarget
target_population = sys.argv[1] # target


print('arg1 =', filename, 'arg2 =', target_population, 'arg3 =', other_population, ' ============')

#folder to which the result will be recorded.
outputfilename = 'outputs'
Constants.TARGED_POPULATION = int(target_population)


Constants.pair_name = '_' + target_population + '_' + other_population
# a list of the two populations numbers (column in fasta) [target, offtarget]
Constants.PAIRWISE_PAIR = [int(e) for e in [target_population, other_population]]

#Runs section descripbed in 'Step 1'
populations_seq_ref = Population_Sequences.read_file_to_sequence_refference(filename)

"""                          -----Step 2-----
    Use 'populations_seq_ref' to find CRISPR/Cas9 30mer target sites in the
    fixed/almost fixed subsequences of a targed population (specied in constants)
    class.

    the thirtymer consists of 4nt + 20nt guide annealing site + 3 nt PAM + 3nt
    the 4nt and 3nt either side of the annealing site are included becouse
    they are used by crisrpon program to calculate annealing effiencey.

    30mers are stored in the biopython sequence format. The .id attribute is the
    scaffold name, and the .description attribute is the index of the the 0th
    element of the subsequence.
"""

thirtymers = Thirtymer.find_30mers(populations_seq_ref)
"""                          -----Step 3-----
    reffine the 30mers found in the target population, keeping only those
    with fixed polymorphims to all other populations.
"""
#for thirtymer in thirtymers:
    #print(thirtymer)
#print('\n\n--------------------')

thirtymers = Thirtymer.refine_to_30mer_with_fixed_dif(thirtymers, populations_seq_ref)
"""                          -----Step 4-----
    construct fasta file for crispron to use
"""
with open(outputfilename + '/30mers' + Constants.pair_name + '.fa', 'w') as f:
    for thirtymer in thirtymers:
        for i in range(len(thirtymer)):
            f.write('>' +  str(thirtymer[i].id) + '_index_' + str(thirtymer[0].index_in_file) + '_population_' + thirtymer[i].population_number +'\n')
            f.write(str(thirtymer[i].seq) + '\n')
    f.close()

with open(outputfilename + '/23mers' + Constants.pair_name + '.fa', 'w') as f:
    for thirtymer in thirtymers:
        for i in range(len(thirtymer)):
            f.write('>' +  str(thirtymer[i].id) + '_index_' + str(thirtymer[0].index_in_file) + '_population_' + thirtymer[i].population_number +'\n')
            f.write(str(thirtymer[i].seq)[4:-3] + '\n')
    f.close()

t1 = time()

"""                          -----Step 4-----
    construct frequenceys file which contains suplmentry data
"""
freq_dict = {'ID': [], 'posstion': [], "polymorphims": [], 'polymorphim_frequency': []}

for thirtymer in thirtymers:
    for i in range(len(thirtymer)):
        if len(thirtymer[i].seq) > 30:
            print(thirtymer[i])
            print('BUG', len(thirtymer[i].seq), '--------------------bug')
            exit()
        freq_dict['ID'].append(str(thirtymer[i].id) + '_index_' + str(thirtymer[0].index_in_file) + '_population_' + thirtymer[i].population_number)
        freq_dict['polymorphims'].append(str(thirtymer[i].polymorphims))
        freq_dict['polymorphim_frequency'].append(str(thirtymer[i].polymorphim_frequency))
        freq_dict['posstion'].append(str(thirtymer[i].reads[0].scaffold_index) + '-' + str(thirtymer[i].reads[-1].scaffold_index))
        

pd.DataFrame(freq_dict).to_csv(outputfilename + '/frequenceys' + Constants.pair_name + '.csv', index = False)




