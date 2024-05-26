import numpy as np

def readfasta(path):
    with open(path, "r") as rf:
        data = rf.read()
        data = data.split(">")[1:]
        formatted_data = list()
        for sequence in data:
            sequence = sequence.split("\n")
            sequence_title = sequence[0]
            sequense_data = "".join(sequence[1:])
            formatted_data.append((sequence_title, sequense_data))
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

    return superlist