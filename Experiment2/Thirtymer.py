from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Population_Sequences import Population_Sequences
from Constants import Constants
from Read import Read
class Thirtymer():
    possistion = ''
    #self.posstion = possition in the file
    def __init__(self, sequence, scaffold, index_in_file, population):
        self.major_allele_freq = 0
        self.polymorphims = []
        self.polymorphim_frequency = []
        self.reads = []
        self.seq = sequence
        self.id = scaffold
        self.index_in_file = index_in_file
        self.population_number = population
        self.complement_to_fasta = False


    def find_30mers( populations_seq_ref):
        thirtymers = []

        for scaffold in populations_seq_ref.keys():
            sequence_record = populations_seq_ref[scaffold].seq_formatted[Constants.PAIRWISE_PAIR.index(Constants.TARGED_POPULATION)] #targeted population at index 0

            if len(sequence_record.seq) >= 30:
                #find PAM site
                for i in range(Constants.BASES_LEFT_OF_PAM, len(sequence_record.seq) - Constants.BASES_RIGHT_OF_PAM):
                    #PAM site found
                    if sequence_record.seq[i] == 'G' and sequence_record.seq[i + 1]  == 'G':
                            # extract 30mer
                            extracted_sequence = sequence_record.seq[i - Constants.BASES_LEFT_OF_PAM: i + Constants.BASES_RIGHT_OF_PAM]
                            # no polymorphic sites included in annealing site and GG on NGG pam site.
                            # these magic numbers refer to the twentymer and  GG of PAM
                            # 4:27 = 23mer + PAM
                            if '-' not in ((extracted_sequence[4:27])[:20] + (extracted_sequence[4:27])[21:]):
                                #save the thirtymer annealing site
                                thirtymer = Thirtymer(extracted_sequence, sequence_record.id, i-Constants.BASES_LEFT_OF_PAM, sequence_record.name)
                                thirtymer.complement_to_fasta = populations_seq_ref[scaffold].complement_to_fasta
                                reads_seq = populations_seq_ref[thirtymer.id].reads[Constants.PAIRWISE_PAIR.index(int(thirtymer.population_number))][i - Constants.BASES_LEFT_OF_PAM: i + Constants.BASES_RIGHT_OF_PAM]
                                # reads object
                                thirtymer.reads = reads_seq
                                # starting possition of the thirtymer.
                                thirtymer.possistion = thirtymer.reads[0].scaffold_index
                                thirtymers.append(thirtymer)
        return thirtymers


    def make_thirymers(current_thirtymer, populations_seq_ref):
        new_thirtymers = []
        for other_population_seq in populations_seq_ref[current_thirtymer.id].seq_formatted:
             if current_thirtymer.population_number != other_population_seq.name:
                 new_seq = other_population_seq.seq[int(current_thirtymer.index_in_file):int(current_thirtymer.index_in_file) + 30]
                 new_thirtymer = Thirtymer(new_seq, other_population_seq.id,
                    current_thirtymer.index_in_file, other_population_seq.name)

                 reads_seq = populations_seq_ref[new_thirtymer.id].reads[Constants.PAIRWISE_PAIR.index(int(new_thirtymer.population_number))][int(current_thirtymer.index_in_file):int(current_thirtymer.index_in_file) + 30]
                 new_thirtymer.reads = reads_seq
                 new_thirtymer.complement_to_fasta = populations_seq_ref[current_thirtymer.id].complement_to_fasta
                 new_thirtymer.possistion = current_thirtymer.possistion
                 #add polymorphim
                 new_thirtymers.append(new_thirtymer)


        for new_thirtymer in new_thirtymers:
            # [4:27] is a bodge fix, as we extract 30mers,however we are only intrested
            # in mismatch of the 23mer.
            for i in range(4,27):
                if new_thirtymer.seq[i] != current_thirtymer.seq[i]:
                    if i not in current_thirtymer.polymorphims:
                        current_thirtymer.polymorphims.append(i)
                        current_thirtymer.polymorphim_frequency.append(current_thirtymer.reads[i].max_base_proportion)
                    new_thirtymer.polymorphims.append(i)
                    new_thirtymer.polymorphim_frequency.append(new_thirtymer.reads[i].max_base_proportion)

        return [current_thirtymer] + new_thirtymers



    # takes 30mers in biopython sequence form and compares them to their
    # equilverlant locus in other populations, if there is fixed difference
    # then they are retained
    def refine_to_30mer_with_fixed_dif(thirtymers, populations_seq_ref):
        refined_30mers = []
        for thirtymer in thirtymers:
            include_30mer = True
            for populations_seq in populations_seq_ref[thirtymer.id].seq_formatted:
                # as above we are only interested in the guide annealing site plus GG of PAM
                offtarget_seq =  populations_seq.seq[int(thirtymer.index_in_file) + 4 :int(thirtymer.index_in_file) + 27]
                offtarget_seq = offtarget_seq[:20] + offtarget_seq[21:]
                # makes sure we dont compare a population to its self.
                if thirtymer.population_number != populations_seq.name:
                    # if consensus sequences are equal.
                    if (thirtymer.seq[4:27])[:20] + (thirtymer.seq[4:27])[21:] == offtarget_seq:
                        include_30mer = False
                if '-' in offtarget_seq:
                    include_30mer = False

                # as sync file is filtered there may be somes gaps
                current_read = thirtymer.reads[0]
                for next_read in thirtymer.reads[1:]:
                    if current_read.scaffold_index - next_read.scaffold_index != 1 and current_read.scaffold_index - next_read.scaffold_index != -1:
                        include_30mer = False
                    current_read = next_read

            if include_30mer:
                refined_30mers.append(Thirtymer.make_thirymers(thirtymer, populations_seq_ref))
        return refined_30mers
