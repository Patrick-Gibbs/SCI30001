import pandas as pd
import numpy as np
import time
import sys

files = [0] # set to range(63) for all populations in dataset
annotation_file = sys.argv[1] # normally 'Lp_Annotation_V1.1.mp.gff3'



scaffold_map = {}
with open(annotation_file, 'r') as f:
    f.readline()
    while 1:
        line = f.readline()
        if not line:
            break

        line = line.split()
        try:
            scaffold_map[line[0]] += [{'type': line[2], 'start': int(line[3]), 'finish': int(line[4])}]
        except KeyError:
            scaffold_map[line[0]] = [{'type': line[2], 'start': int(line[3]), 'finish': int(line[4])}]

hits = []
for i in files:
    print(i)
    with open('Results/Results_' + str(i) + '.csv', 'r') as f:
        f.readline()
        while 1:
            line = f.readline()
            if not line:
                break
            line = line.split(',')
            print(line[1])
            hash = line[1][:line[1].index('-')]

            start_finish = sorted(line[2].split('-'))
            start = int(start_finish[0])
            finish = int(start_finish[1])

            no_hit = 1
            try:
                isInGene = False
                for e in scaffold_map[hash]:
                    if finish > e['start'] and start < e['finish']:
                        if e['type'] == 'gene':
                            isInGene = True
                hits.append(isInGene)
            except KeyError:
                hits.append(False)

    data = pd.read_csv('Results/Results_' + str(i) + '.csv')
    data.insert(6, "in_gene", hits, allow_duplicates = False)
    data.to_csv('Results' + str(i) + '.csv', index=False)
    hits = []
