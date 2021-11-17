from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Population_Sequences import Population_Sequences
from Constants import Constants
from Read import Read

import copy
class AnnealingSite():
    possistion = ''
    #self.posstion = possition in the file
    def __init__(self, sequence, scaffold, index_in_file, population):
        # frequencey of the major allel in the an annealing site
        self.reads = []
        self.seq = sequence
        self.id = scaffold
        self.index_in_file = index_in_file
        self.population_number = population
        self.complement_to_fasta = False

    """ finds all annealing sites where differentation between the targeted
        and offtarget populations is above a set threshold

        inputs: populations_seq_ref, a dictionary containing two keys, the +
        and - strands of scaffolds scaffolds. values few each key is infomation
        of the sequence.

        outputs: annealing_sites, valid annealing sites where differentation
        between the target and offtarget are below a set threshold.
    """
    def find_annealing_sites(populations_seq_ref):
        # for each scaffold valid annealing sites are stored in this list
        annealing_sites = []

        # scaffold name for possitive and negative strand.
        for scaffold in populations_seq_ref.keys():
            sequence_record = populations_seq_ref[scaffold].seq_formatted[Constants.TARGED_POPULATION]
            if len(sequence_record.seq) >= 23: # 20mer guide + PAM
                #find PAM site
                for i in range(Constants.BASES_LEFT_OF_PAM,len(sequence_record.seq) - Constants.BASES_RIGHT_OF_PAM):
                    # the first element in the PAM site (index 20) does not matter.
                    """if i == 20: #population 20 was found to be an outlyer
                        continue"""
                    #PAM site found
                    if sequence_record.seq[i] == 'G' and sequence_record.seq[i + 1]  == 'G':
                        include_this_pam = True
                        # extranct the consensus sequences of the target
                        target_annealing_site = (sequence_record.seq[i - Constants.BASES_LEFT_OF_PAM: i
                            + Constants.BASES_RIGHT_OF_PAM + 1])
                        # extract exact reads for each possition in the target
                        extracted_reads_target = (populations_seq_ref[scaffold].reads[Constants.TARGED_POPULATION][
                            i - Constants.BASES_LEFT_OF_PAM: i + Constants.BASES_RIGHT_OF_PAM + 1])

                        if '-' in target_annealing_site[:20] + target_annealing_site[:21 ]: break # we dont include if it contains bad reads.

                        # ensure no gaps, this means the program can read filtered sync files
                        current_read = extracted_reads_target[0]
                        for next_read in extracted_reads_target[1:]:
                            if current_read.scaffold_index - next_read.scaffold_index != 1 and current_read.scaffold_index - next_read.scaffold_index != -1:
                                include_this_pam = False
                                break
                            current_read = next_read
                        if not include_this_pam: #this makes it alittle faster
                            break

                        # extract all reads in the targeted annealing site for every population
                        extracted_reads = [reads[i - Constants.BASES_LEFT_OF_PAM: i + Constants.BASES_RIGHT_OF_PAM + 1] for reads in populations_seq_ref[scaffold].reads]

                        # gets the frequencey for any single nucleotide in a population
                        def freq_of_NT(nt, read):
                            if nt == 'A': freq = read.a_reads
                            if nt == 'T': freq = read.t_reads
                            if nt == 'C': freq = read.c_reads
                            if nt == 'G': freq = read.g_reads
                            return freq/read.total_reads

                        # for every nucelotide possition in the extracted annealing site
                        for read_index in range(len(target_annealing_site)):
                            # stores the frequencey of each targeted nucelotide in each population
                            frequenceys = []
                            # as above by the key serves as a refference to the population the frequencey came.
                            freq_dict = {}
                            # read object for the target population at the current possition.
                            target_reads = extracted_reads[Constants.TARGED_POPULATION][read_index]
                            # consensus basepair for the target population at this possition.
                            max_freq_target = target_reads.max_base_proportion
                            target_nt = target_reads.max_base

                            # iterate every population other population
                            for pop_num in range(Constants.POPULATION_NUMBER):
                                # 20 is and outlyer.
                                if pop_num != Constants.TARGED_POPULATION:
                                    # check reads in every population are of good quality.
                                    if '-' == extracted_reads[pop_num][read_index].max_base:
                                        possition_score = 0
                                        include_this_pam = False
                                        break

                                    # records the frequencey of each nt.
                                    frequenceys.append(freq_of_NT(target_nt, extracted_reads[pop_num][read_index]))
                                    freq_dict[pop_num] = freq_of_NT(target_nt, extracted_reads[pop_num][read_index])

                            # average difference
                            if include_this_pam: possition_score = max_freq_target - sum(frequenceys)/len(frequenceys)
                            if possition_score > Constants.MIN_DIFFERENCE_SCORE:
                                # make the object
                                # target_sequence, scaffold, index_in_file, population
                                annealing_site = AnnealingSite(target_annealing_site, sequence_record.id,
                                i-Constants.BASES_LEFT_OF_PAM, sequence_record.name)
                                annealing_site.polymorphim_index  = read_index
                                annealing_site.complement_to_fasta = populations_seq_ref[scaffold].complement_to_fasta
                                annealing_site.reads = extracted_reads[Constants.TARGED_POPULATION]
                                annealing_site.possistion = annealing_site.reads[read_index].scaffold_index
                                annealing_site.max_freq_target = max_freq_target
                                annealing_site.possition_score = possition_score
                                annealing_site.frequenceys = frequenceys
                                annealing_site.freq_dict = copy.deepcopy(freq_dict)
                                annealing_sites.append(annealing_site)
        return annealing_sites
