from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
class Read:
    """ this object holds all infomation for a locus specific to one population"""
    def __init__(self, scaffold, scaffold_index, sync_infomation, population_number):
        #infomation about secquence

        sync_infomation_list = [int(x) for x in sync_infomation.split(':')]
        # number of reads for each NC at a locus for a population
        self.a_reads = sync_infomation_list[0]
        self.t_reads = sync_infomation_list[1]
        self.c_reads = sync_infomation_list[2]
        self.g_reads = sync_infomation_list[3]
        self.total_reads = sum(sync_infomation_list)


        #scaffold and scaffold_index
        self.scaffold = scaffold


        self.population_number = population_number

        self.scaffold_index = int(scaffold_index)


        if self.total_reads < 50:
            self.max_base = '-'
            self.max_base_proportion = 0
            return
        # proportions of reads
        self.base_proportions = self.get_proportions()
        self.max_base_proportion = max(self.base_proportions)

        #index of the most common base pair in the popuation
        index = self.base_proportions.index(self.max_base_proportion)
        if index == 0:
            self.max_base = 'A'
        elif index == 1:
            self.max_base = 'T'
        elif index == 2:
            self.max_base = 'C'
        elif index == 3:
            self.max_base = 'G'

    def complement(self):
        temp_a = self.a_reads
        temp_t = self.t_reads
        temp_c = self.c_reads
        temp_g = self.g_reads

        self.a_reads = temp_t
        self.t_reads = temp_a
        self.c_reads = temp_g
        self.g_reads = temp_c
        return self

    # used to get the proportions of at an index for a population
    def get_proportions(self):
        proportions = []
        try:
            proportions.append(self.a_reads/self.total_reads)
            proportions.append(self.t_reads/self.total_reads)
            proportions.append(self.c_reads/self.total_reads)
            proportions.append(self.g_reads/self.total_reads)
        except:
            print("BUG reads =====", self.scaffold, self.population_number)
            exit()
        return proportions
