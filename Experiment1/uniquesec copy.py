from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

from Population_Sequences import Population_Sequences
from Constants import Constants
from Read import Read
from time import time

import sys



t0 = time()
filename = sys.argv[1]
target_population = sys.argv[2]
outputfilename = 'outputs'
Constants.TARGED_POPULATION = int(target_population)
annealing_sites = Population_Sequences.read_file_to_sequence_refference(filename)

results_dict = {"scafoldID": [], "possition": [], "exact_possition": [], "polymorphim_index": [],  'target_freq': [], "spec_score": [], 'seq':[], "target_population": []}
for i in range(Constants.POPULATION_NUMBER):
    if i != Constants.TARGED_POPULATION:
        results_dict[i] = []

for annealing_site in annealing_sites:
    results_dict["scafoldID"].append(annealing_site.id)
    results_dict["exact_possition"].append(annealing_site.possistion)
    results_dict["possition"].append(str(annealing_site.reads[0].scaffold_index) + '-' + str(annealing_site.reads[-1].scaffold_index))
    results_dict["polymorphim_index"].append(annealing_site.polymorphim_index)
    results_dict['target_population'].append(Constants.TARGED_POPULATION)
    results_dict["spec_score"].append(annealing_site.possition_score)
    results_dict['seq'].append(str(annealing_site.seq))
    results_dict["target_freq"].append(annealing_site.max_freq_target)

    for i in range(Constants.POPULATION_NUMBER):
        if i != Constants.TARGED_POPULATION:
            results_dict[i].append(annealing_site.freq_dict[i])

results = pd.DataFrame(results_dict)
results.to_csv("Results_" + str(Constants.TARGED_POPULATION) + '.csv')
print(results)
