from Constants import Constants
from Read import Read
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

    def add_to_sequences(self, input_line):
        """ as each line of the sync file is read the methord adds one nucleotide
            to each current population sequence """
            # reads in just the target and off target sequence
        for i in [j + Constants.SEQUENCE_INFO_OFFSET for j in Constants.PAIRWISE_PAIR]:
            read = Read(input_line[0], input_line[1], input_line[i], i - Constants.SEQUENCE_INFO_OFFSET)
            # here we are building a sequnce one read at a time, appending to two lists only.
            # .index is used as it is a safer than using 0,1 bc target and off target wont get confused.
            self.reads[Constants.PAIRWISE_PAIR.index(i - Constants.SEQUENCE_INFO_OFFSET)].append(read)
            # checks that each read is not greater than a set minium threshold for a considered consusus
            # currently this is set to 0, so the max freq will always be valid and considered the consusus.
            if read.max_base_proportion > Constants.MAX_BASE_PROPORTION:
                self.sequences[Constants.PAIRWISE_PAIR.index(i - Constants.SEQUENCE_INFO_OFFSET)].append(read.max_base)
            else:
                self.sequences[Constants.PAIRWISE_PAIR.index(i - Constants.SEQUENCE_INFO_OFFSET)].append('-')

    """ converts all sequences currently stored as list in this object
    to biopython seq format. rational is sequences are built as lists
    and then stored as arrays """
    def to_biopython_seq_format(self, current_scaffold):
        #converts to biopython sequence format
        for i in range(len(self.sequences)):
            #converts list of chars to biopython sequnce format
            self.seq_formatted.append(SeqRecord(Seq("".join(self.sequences[i])),
                id=current_scaffold, name=str(Constants.PAIRWISE_PAIR[i])))

    def to_biopython_seq_format_reversed(self, current_scaffold):
        #converts to biopython sequence format
        for i in range(len(self.sequences)):
            # converts list of chars to biopython sequnce format reverses order and
            # changes to complement
            #+ '-s' + str(Population_Sequences.fasta_index_for_zeroth)
            self.seq_formatted.append(SeqRecord(Seq("".join(self.sequences[i][::-1])).complement(),
                id=current_scaffold, name=str(Constants.PAIRWISE_PAIR[i])))
            # order is reversed for reads, included .complment method can also be used.
            # however it is not done for here for speed speed as nucelotides are not used from
            # reads later in the program.
            self.reads[i] = self.reads[i][::-1]

    """
    adds biopython version of the sequence to seq ref, including the reverses
    complement.
    """
    def update_seq_ref(populations_seq_ref, population_sequences, current_scaffold):
        population_sequences.to_biopython_seq_format(current_scaffold + "-p")
        populations_seq_ref[current_scaffold + "-p"] = population_sequences
        population_sequences = copy.deepcopy(population_sequences)
        population_sequences.seq_formatted = []
        population_sequences.complement_to_fasta = True
        population_sequences.to_biopython_seq_format_reversed(current_scaffold + "-m")
        populations_seq_ref[current_scaffold + "-m"] = population_sequences
        return populations_seq_ref

    def read_file_to_sequence_refference(filename):
        """ reads sync file and returns a refference encaputing the consensus
            sequence for each population under each refference """

        #read in sync file
        file = open(filename, 'r')
        line = file.readline()
        line = line.split()
        current_scaffold = line[0]

        # for each refference the consensus of all population is kept track
        # of in population sequences object
        population_sequences = Population_Sequences(line)

        # dictionary where scafold name is key and items are lists of sequence
        # for each population
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
                #print(population_sequences.sequences)
                populations_seq_ref = Population_Sequences.update_seq_ref(populations_seq_ref, population_sequences, current_scaffold)
                current_scaffold = line[0]
                population_sequences = Population_Sequences(line)
                population_sequences.add_to_sequences(line)
            line = file.readline()
        populations_seq_ref = Population_Sequences.update_seq_ref(populations_seq_ref, population_sequences, current_scaffold)

        file.close()
        return populations_seq_ref
