import sys
import numpy as np

'''This is the main mehtod that is the diver of the program'''
def main():
    seq= None
    score = None
    penalty = None
    # Checks for valid inputs 
    try:
        seq = sys.argv[1] 
        score = sys.argv[2]
        penalty = int(sys.argv[3])
    except: 
        print("An error as occured when reading in the inputs")
        return

    # obtain the sequences and scoring matrix from the files
    s_matrix = None
    all_seq = None
    with open(seq) as file:
        lines = file.readlines()
        all_seq = [line.rstrip().upper() for line in lines]
        file.close()
    with open(score) as file:
        lines = file.readlines()
        s_matrix = [line.rstrip().split() for line in lines]
        file.close()

    # Check for invalid inputs and print out error statements 
    if len(all_seq) != 2:
        print("ERROR: The input has the incorrect number of sequences to align")
        return
    valid_char_1 = s_matrix[0][1:]
    valid_char_2 = list(np.array(s_matrix)[:, 0])[1:]
    for i in valid_char_1:
        if valid_char_1.count(i) > 1:
            print("ERROR: Your matrix has a double count")
        if i not in  valid_char_2:
            print("ERROR: Your scoring scheme does not match")
            return
    for i in valid_char_2:
        if valid_char_2.count(i) > 1:
            print("ERROR: Your matrix has a double count")
        if i not in  valid_char_1:
            print("ERROR: Your scoring scheme does not match")
            return
    for i in all_seq[0]:
        if i not in valid_char_1:
            print("ERROR: You have an invalid charater")
            return
    for i in all_seq[1]:
         if i not in valid_char_1:
            print("ERROR: You have an invalid charater")
            return

    # Finds the optimal alignment
    optimal_aligment(all_seq, s_matrix, penalty)

'''This method finds the optimal aligment between two sequences'''
def optimal_aligment(all_seq, s_matrix, penalty):
    seq_1 = all_seq[0]
    seq_2 = all_seq[1]
    matrix = []
    trace = []

    # Create db table as well as backtraking matrix  
    for i in range(len(seq_1)+1):
        m_subarray = []
        t_subarray = []
        for j in range(len(seq_2)+1):
            m_subarray.append(0)
            t_subarray.append('0')
        matrix.append(m_subarray)
        trace.append(t_subarray)

    # set up matrix borders of the DP table
    for i in range(1, len(seq_1)+1):
        matrix[i][0] = i*penalty
    for i in range(1, len(seq_2)+1):
        matrix[0][i] = i*penalty

    # set up matrix borders of the backtrace table
    for i in range(1, len(seq_1)+1):
        trace[i][0] = " up "
    for i in range(1, len(seq_2)+1):
        trace[0][i] = "left"
    trace[0][0] = "done"

    # run through the entire matrix and fill out the DP table and backtracking 
    for i in range(1, len(seq_1)+1):
        for j in range(1, len(seq_2)+1):
            left = matrix[i][j-1] + penalty
            up = matrix[i-1][j] + penalty
            diag = matrix[i-1][j-1] + find_score(s_matrix, seq_1[i-1], seq_2[j-1])
            
            matrix[i][j] = max(left, up, diag)
            
            if matrix[i][j] == left:
                trace[i][j] = "left"
            if matrix[i][j] == up:
                trace[i][j] = " up "
            if matrix[i][j] == diag:
                trace[i][j] = "diag"

    # get the sequence alignment 
    seq_1, seq_2 = getalign(seq_1, seq_2, trace)
    print(seq_1)
    print(seq_2)
    print(matrix[-1][-1])

'''This method finds the match and mismatch score'''
def find_score(s_matrix, a, b):
    a_index = s_matrix[0].index(a)
    b_index = list(np.array(s_matrix)[:,0]).index(b)
    output = int(s_matrix[a_index][b_index])
    return output

'''This method gets the alignment sequences of both'''
def getalign(seq_1, seq_2, trace):
    x_seq = ""
    y_seq = ""
    x_len = len(seq_1)
    y_len = len(seq_2)
    while (x_len > 0 or y_len > 0):
        if trace[x_len][y_len] == "diag":
            x_seq = (seq_1[x_len-1]) + x_seq
            y_seq = (seq_2[y_len-1]) + y_seq
            x_len -= 1
            y_len -= 1
        elif trace[x_len][y_len] == "left":
            x_seq = ("-") + x_seq
            y_seq = (seq_2[y_len-1]) + y_seq
            y_len -= 1
        elif trace[x_len][y_len] == " up ":
            x_seq = (seq_1[x_len-1]) + x_seq
            y_seq = ("-") + y_seq
            x_len -= 1
        elif trace[x_len][y_len] == "done":
            break
    return x_seq, y_seq

if __name__ == "__main__":
    main()