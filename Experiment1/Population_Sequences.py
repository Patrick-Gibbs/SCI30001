from Constants import Constants
from Read import Read
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from AnnealingSite import AnnealingSite

import copy

class Population_Sequences():
    def __init__(self, first_line):
         """ stores the sequences for each population for refference
         format: [[population 1 sequnce], [population 2 sequence], ...] """
         self.fasta_index_for_zeroth = int(first_line[1])
         self.seq_formatted = []
         self.sequences = [[] for x in range(Constants.POPULATION_NUMBER)]
         self.reads = [[] for x in range(Constants.POPULATION_NUMBER)]
         self.complement_to_fasta = False




    """ as each line of the sync file is read the methord adds one nucleotide
        to each current population sequence """
    def add_to_sequences(self, input_line):
        for i in [j + Constants.SEQUENCE_INFO_OFFSET for j in range(Constants.POPULATION_NUMBER)]:
            read = Read(input_line[0], input_line[1], input_line[i], i - Constants.SEQUENCE_INFO_OFFSET)

            self.reads[i - Constants.SEQUENCE_INFO_OFFSET].append(read)
            if read.max_base_proportion > Constants.MAX_BASE_PROPORTION:
                self.sequences[i - Constants.SEQUENCE_INFO_OFFSET].append(read.max_base)
            else:
                self.sequences[i - Constants.SEQUENCE_INFO_OFFSET].append('-')

    """ converts all sequences currently stored as list in this object
        to biopython seq format. rational is sequences are built as lists
        and then stored as arrays – arrays will iterate faster. """
    def to_biopython_seq_format(self, current_scaffold):

        #converts to biopython sequence format
        for i in range(len(self.sequences)):
            #converts list of chars to biopython sequnce format
            self.seq_formatted.append(SeqRecord(Seq("".join(self.sequences[i])),
                id=current_scaffold, name=str(i)))

    """ Same as the above method however stores the reverse complement"""
    def to_biopython_seq_format_reversed(self, current_scaffold):
        #converts to biopython sequence format
        for i in range(len(self.sequences)):
            # converts list of chars to biopython sequnce format reverses order and
            # changes to complement
            self.seq_formatted.append(SeqRecord(Seq("".join(self.sequences[i][::-1])).complement(),
                id=current_scaffold, name=str(i)))
            self.reads[i] = self.reads[i][::-1]



    def update_seq_ref(populations_seq_ref, population_sequences, current_scaffold):
        population_sequences.to_biopython_seq_format(current_scaffold + "-p")
        populations_seq_ref[current_scaffold + "-p"] = population_sequences
        population_sequences = copy.deepcopy(population_sequences)
        # this is an attempt to save some memory.
        population_sequences.seq_formatted = []
        population_sequences.complement_to_fasta = True
        population_sequences.to_biopython_seq_format_reversed(current_scaffold + "-m")
        populations_seq_ref[current_scaffold + "-m"] = population_sequences
        return populations_seq_ref

    def read_file_to_sequence_refference(filename):
        """ reads sync file and returns a refference encaputing the consensus
            sequence for each population under each refference """


        annealing_site_objects = []

        #read in sync file
        file = open(filename, 'r')
        line = file.readline()
        line = line.split()
        current_scaffold = line[0]

        # by defult all popoualtions are incldued
        if not Constants.POPULATION_NUMBER:
            Constants.POPULATION_NUMBER = len(line) - 3 # minus 3 bc first 3 rows are not pops

        # for each refference the consensus of all population is kept track
        # of in population sequences object
        population_sequences = Population_Sequences(line)

        # dictionary where scafold name is key and items are lists of sequences
        # for each population.
        populations_seq_ref = {}
        # read sync line by line
        while(line):
            # line is already split on first iteration
            if type(line) != type([]):
                line = line.split()
            # we keep track of what scaffold we are in
            if line[0] == current_scaffold:
                population_sequences.add_to_sequences(line)
            # new scaffold
            else:
                populations_seq_ref = Population_Sequences.update_seq_ref(populations_seq_ref, population_sequences, current_scaffold)
                # find all anealing sites in the scafold we just read in.
                annealing_site_objects += AnnealingSite.find_annealing_sites(populations_seq_ref)
                populations_seq_ref = {}
                current_scaffold = line[0]
                population_sequences = Population_Sequences(line)
                population_sequences.add_to_sequences(line)
            line = file.readline()
        populations_seq_ref = Population_Sequences.update_seq_ref(populations_seq_ref, population_sequences, current_scaffold)

        file.close()
        return annealing_site_objects
