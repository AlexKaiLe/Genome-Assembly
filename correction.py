import sys
import numpy as np
# from progress.bar import Bar

'''Obtains the frequent kmers from all the reads'''
def gen_freq_kmers(reads, k, t):
    kmer_freq = {}
    freq_kmers = []
    # run through all the reads
    for i in range(len(reads)):
        read = reads[i]
        # run through the length of a read
        for j in range(len(read)-k+1):
            cur_kmer = read[j:j+k]
            # add values to dict or increment values of dict
            if cur_kmer not in kmer_freq:
                kmer_freq[cur_kmer] = 1
            else:
                kmer_freq[cur_kmer] += 1
    # run through the dict and if the freq is above or equal to the thresh, add it to the frequent k-mer list
    for kmer in kmer_freq:
        if kmer_freq[kmer] >= t:
            # print(t)
            freq_kmers.append([kmer_freq[kmer], kmer])
    return kmer_freq, sorted(freq_kmers, reverse= True)

'''Obtains the infrequent kmers and the index for a read'''
def get_infreq(read_kmer, freq_list):
    infreq = []
    infreq_index = []
    # run through the kmers of the read
    for i in range(len(read_kmer)):
        # if the kmer is not in the freq list, than add it to a infreq list
        if read_kmer[i] not in freq_list:
            infreq.append(read_kmer[i])
            infreq_index.append(i)
    return infreq, infreq_index

'''Obtains kmers from a sequence'''
def get_kmers(seq, k):
    kmers = []
    index = []
    # run through a sequence to generate kmers
    for i in range(len(seq)-k+1):
        kmers.append(seq[i:i+k])
        index.append(index)
    return kmers, index

'''Obtain a list of freq kmers that best replaces the infreq kmer'''
def optimal_replace(infreq_kmer, freq_kmers, d, k):
    replace_kmer = []
    # run through all the freq kmers
    for fq_kmer in freq_kmers:
        diff = 0
        # run through the length of k
        for i in range(k):
            # add to difference score if the freq kmer is not the same
            if infreq_kmer[i] != fq_kmer[i]:
                diff += 1
        # if the hamming distnace is lower or equal to d, then add it to infrequent kmers
        if diff <= d:
            replace_kmer.append([diff, fq_kmer])
    return sorted(replace_kmer, key=lambda x: x[1])

'''Correct all the incorrect reads'''
def correct(reads, freq_kmers, d, k):
    freq_kmers = np.array(freq_kmers)[:,1]
    corrected_reads = []
    corrected_index = []
    replacable = {}

    # run through all the reads in the file
    # for i in Bar('Processing').iter(range(len(reads))):
    for i in range(len(reads)):
        # print(i/len(reads))
        cur_read = reads[i]
        loop_count = 0

        # keep on altering the read until it reaches a limit or we are satisfied 
        while True:
            # get the kmers of the current read and the infreq kmers that appear
            read_kmer, kmer_index = get_kmers(cur_read, k)
            infreq, infreq_index = get_infreq(read_kmer, freq_kmers)
            # If there are infreq kmers, replace them
            if len(infreq) > 0:
                new_list = []
                # run through the list of infreq kmers
                for j in range(len(infreq)):
                    replace_kmer = infreq[j]
                    replace_index = infreq_index[j]
                    # find a replacable value for the kmer of interest
                    if replace_kmer not in replacable:
                        new_kmer = optimal_replace(replace_kmer, freq_kmers, d, k)
                        if len(new_kmer) == 0:
                            continue
                        replacable[replace_kmer] = new_kmer[0][1]
                    new_kmer = replacable[replace_kmer]
                    # create new replaceable sequence
                    new_read = cur_read[:replace_index] + new_kmer + cur_read[replace_index+k:]

                    # check how many k-mers this new sequence introduces
                    new_read_kmer, new_kmer_index = get_kmers(new_read, k)
                    new_infreq, new_infreq_index = get_infreq(new_read_kmer, freq_kmers)
                    new_list.append([len(new_infreq), replace_index, new_read])

                # if there are corrected sequences, choose the one that creates the least amount of k-mers
                if len(new_list) > 0:
                    cur_read = sorted(new_list)[0][2]
                # add the index of the changed read to a list
                if str(i) not in corrected_index:
                    corrected_index.append(str(i))
            else:
                break

            if loop_count > 20:
                break
            loop_count += 1
        
        corrected_reads.append(cur_read)
    return corrected_reads, corrected_index

'''method only used for the application'''
def HIV(reads, k, t, d):
    print(k, t, d)
    kmer_freq, freq_kmers = gen_freq_kmers(reads, k, t)
    print("Number of frequent k-mers:", len(freq_kmers), "/", len(kmer_freq), "=", len(freq_kmers)/len(kmer_freq))
    print(freq_kmers[0], freq_kmers[-1])
    corrected_reads, corrected_index = correct(reads, freq_kmers, d, k)
    print("Number of correted reads", len(corrected_reads))
    print('--------------------')
    return corrected_reads