import numpy as np

def AA(codon):
    codon_table = {"UUU":"F",
               "UUC":"F",
               "UUA":"L",
               "UUG":"L",
               "UCU":"S",
               "UCC":"S",
               "UCA":"S",
               "UCG":"S",
               "UAU":"Y",
               "UAC":"Y",
               "UAA":"stop",
               "UAG":"stop",
               "UGU":"C",
               "UGC":"C",
               "UGA":"stop",
               "UGG":"W",
               "CUU":"L",
               "CUC":"L",
               "CUA":"L",
               "CUG":"L",
               "CCU":"P",
               "CCC":"P",
               "CCA":"P",
               "CCG":"P",
               "CAU":"H",
               "CAC":"H",
               "CAA":"Q",
               "CAG":"Q",
               "CGU":"R",
               "CGC":"R",
               "CGA":"R",
               "CGG":"R",
               "AUU":"I",
               "AUC":"I",
               "AUA":"I",
               "AUG":"M",
               "ACU":"T",
               "ACC":"T",
               "ACA":"T",
               "ACG":"T",
               "AAU":"N",
               "AAC":"N",
               "AAA":"K",
               "AAG":"K",
               "AGU":"S",
               "AGC":"S",
               "AGA":"R",
               "AGG":"R",
               "GUU":"V",
               "GUC":"V",
               "GUA":"V",
               "GUG":"V",
               "GCU":"A",
               "GCC":"A",
               "GCA":"A",
               "GCG":"A",
               "GAU":"D",
               "GAC":"D",
               "GAA":"E",
               "GAG":"E",
               "GGU":"G",
               "GGC":"G",
               "GGA":"G",
               "GGG":"G"}
    return codon_table[codon]


def Translation(rna):
    protein = np.empty((1, 0))
    for i in range(0, len(rna), 3):
        protein = np.append(protein, AA(rna[i:i+3]))

    return "".join(protein.tolist())


def reverse_complement(dna):
    complement_map = {"G":"C", "A":"T", "T":"A","C":"G"}
    complement = np.array([complement_map[b] for b in dna])
    return "".join(complement[::-1].tolist())

def Transkription(dna):
    return reverse_complement(dna).replace("T", "U")


class FASTA:
    def __init__(self, tag, sequence):
        self.tag = tag
        self.sequence = sequence
        self.rev_complement = reverse_complement(sequence)
        self.RNA = sequence.replace("T", "U")
        self.RNA_rev_complement = reverse_complement(sequence).replace("T", "U")


def readfasta(path):
    with open(path, "r") as rf:
        data = rf.read()
        data = data.split(">")[1:]
        formatted_data = np.empty((1, 0))
        for sequence in data:
            sequence = sequence.split("\n")
            sequence_title = sequence[0]
            sequence_data = "".join(sequence[1:])
            formatted_data = np.append(formatted_data, [FASTA(sequence_title, sequence_data)]) # type: ignore
        return formatted_data


def align_reads(data_list):
    #from LONG exercise
    seq_list = data_list
    #superlist is an array with the col0 as the labels, col1 is the corresponding seq, 
    #col2 index of the end of the left fragment

    superlist = np.empty((0, 3))


    #format data in halfes to check for overlaps:
    #col0 is label, col1 is seq
    seq_list = np.array([[n[0], n[1]] for n in seq_list])
    # print(len(seq_list))


    superlist = np.append(superlist, [[seq_list[0, 0], seq_list[0, 1], 0]], axis=0)
    print((superlist))
    seq_list = seq_list[1::, ::]

    change = True
    while change:

        change = False

        right_end = superlist[-1, 1][len(superlist[-1, 1])//2:]

        for row in range(0, len(seq_list)):
            if right_end in seq_list[row, 1]:
                label = seq_list[row, 0]
                seq = seq_list[row, 1]
                for n in range(len(right_end), len(seq)):
                    if right_end == seq[n-len(right_end):n]:
                        superlist = np.append(superlist, [[label, seq, n]], axis=0)
                        seq_list = np.delete(seq_list, row, axis=0)
                        change = True
                if change:
                    break
                else:
                    print("Alignment found. Couldn't be processed!")

    #left end:

    change = True
    while change:

        change = False

        left_end = superlist[0, 1][:len(superlist[-1, 1])//2]

        for row in range(0, len(seq_list)):
            if left_end in seq_list[row, 1]:
                label = seq_list[row, 0]
                seq = seq_list[row, 1]
                for n in range(len(left_end), len(seq)):
                    if left_end == seq[-n: -n+len(left_end)]:
                        superlist = np.append([[label, seq, 0]], superlist, axis=0)
                        superlist[1, 2] = n
                        seq_list = np.delete(seq_list, row, axis=0)
                        change = True
                if change:
                    break
                else:
                    print("Alignment found. Couldn't be processed!")

    #actually align
    superstring = superlist[0, 1]
    for read in superlist[1:, ]:
        superstring += read[1][int(read[2]):]

    return superstring

def find_ORF(data):
    #from ORF Problem
    stop_codons = ["UAA", "UAG", "UGA"]

    reading_frames = np.empty((1, 0))

    for seq in [data.RNA, data.RNA_rev_complement]:

        for n in range(0, len(seq)):

            if seq[n:n+3] == "AUG":
                for i in range(n+3, len(seq), 3):

                    if seq[i:i+3] in stop_codons:
                        reading_frames = np.append(reading_frames, seq[n:i])
                        break

    return reading_frames


def find_motifs(seq, subseq):    

    loc = np.empty((1, 0))

    for x in range(0, len(seq)):

        if seq[x:x+len(subseq)] == subseq:
            loc = np.append(loc, x)

    return loc
