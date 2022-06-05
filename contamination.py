import sys

'''Obtain the inputs from the sh script'''
def obtain_inputs():
    try:
        # Checks for valid inputs 
        r = sys.argv[1] 
        v = sys.argv[2]
        k = sys.argv[3]

        # obtain the sequences and scoring matrix from the files
        with open(r) as file:
            lines = file.readlines()
            reads = [line.rstrip().upper() for line in lines]
            file.close()

        with open(v) as file:
            vector = file.readlines()[0].rstrip().upper()
            file.close()
        
        k = int(k)
        return reads, vector, k

    except:
        return None, None, None

'''Create kmers given a sequence'''
def create_kmers(vector, k):
    k_mers = []
    for i in range(len(vector)-k+1):
        k_mers.append([i, vector[i:i+k]])
    return k_mers

'''Run through all the reads and remove the contamination'''
def remove_contam(reads, k_mers, vector, k):
    clean = []
    indices = []
    # Run through all the reads
    for i in range(len(reads)):
        read = reads[i]
        far_left = []
        far_right = []
        # run through all the k_mers of the vector
        for j in range(len(k_mers)):
            index = k_mers[j][0]
            k_mer = k_mers[j][1]
            # obtain the left and right most indicies to remove 
            far_left.append(front(read, k_mer, index, vector, k))
            far_right.append(back(read, k_mer, index, vector, k))
        
        #take the right most left index 
        max_left = max(far_left)
        #take the left most right index
        min_right = min(far_right)
        # check if the bounds extend the range of the read 
        if max_left > 0 or min_right < len(read):
            indices.append(i)

        # append all cleaned reads to a list
        clean.append(read[max_left:min_right])
    return indices, clean

'''Get the contamination at the front'''
def front(read, k_mer, index, vector, k):
    left = 0
    # forward
    if check_match(read[0:k], k_mer):
        # run through the read and extend from the left to the right
        for i in range(len(read)):
            vector_index = index + i
            read_index = i

            # check if the value of the read and vector are the same 
            if vector_index >= len(vector):
                break
            if read[read_index] == vector[vector_index]:
                left = read_index + 1
            else:
                break
    return left

'''Get the contamination at the end'''
def back(read, k_mer, index, vector, k):
    right = len(read)
    # backward
    if check_match(read[-k:-1], k_mer):
        # run through the read and extend from the right to the left
        for i in range(len(read)):
            vector_index = index - i
            read_index = len(read)-k-i

            # check if the value of the read and vector are the same 
            if vector_index < 0:
                break
            if read[read_index] == vector[vector_index]:
                right = read_index
            else:
                break
    return right

'''Check if the prefix/sufix are the same'''
def check_match(read, k_mer):
    if len(read) != len(read):
        return False
    for i in range(len(read)):
        if k_mer[i] != read[i]:
            return False
    return True

'''Print the output of the contamination'''
def output(indices, clean):
    print(",".join(indices))
    print("--------------------")
    for i in range(len(clean)):
        print(clean[i])

def contamination(reads, vector, k):
    # Check if there are any errors with reading in the data
    indices, clean = "", reads
    try:
        if k < len(reads[0]) and k < len(vector):
            k_mers = create_kmers(vector, k)
            indices, clean = remove_contam(reads, k_mers, vector, k)
        # output(indices, clean)
        return indices
    except:
        # output(indices, clean)
        return indices
