import sys
import contamination
import correction
import debruijn

'''Get inputs from sh file'''
def obtain_inputs():
    # Checks for valid inputs 
    reads = sys.argv[1]
    vector = sys.argv[2]

    contam_k = int(sys.argv[3])
    
    kmer_length = int(sys.argv[4])

    correct_k = int(sys.argv[5])
    correct_t = int(sys.argv[6])
    correct_d = int(sys.argv[7])

    with open(reads) as file:
        lines = file.readlines()
        reads = [line.rstrip().upper() for line in lines]
        file.close()

    with open(vector) as file:
        vector = file.readlines()[0].rstrip().upper()
        file.close()
    
    
    return reads, vector, [contam_k, kmer_length, correct_k, correct_t, correct_d]

'''From the vector, create kmers'''
def create_kmers(vector, k, k_mers):
    for i in range(len(vector)-k+1):
        k_mers.append(vector[i:i+k])

'''Preform contamination, correction and debruijn'''
def application(reads, vector, numb):
    # contamination
    contam_indices = contamination.contamination(reads, vector, numb[0])
    print("Number of contaminated reads:", len(contam_indices))
    print("Number of uncontaminated:", len(reads) - len(contam_indices))

    for i in reversed(sorted(contam_indices)):
        reads.pop(i)

    # correction
    corrected_reads = correction.HIV(reads, numb[2], numb[3], numb[4])

    kmers = []
    for i in range(len(corrected_reads)):
        create_kmers(reads[i], numb[1], kmers)
        
    # debruijn
    return debruijn.application(kmers)

def main():
    reads, vector, numb = obtain_inputs()
    application(reads, vector, numb)

if __name__ == "__main__":
    main()