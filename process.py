from Bio import SeqIO
import numpy as np
import sys
import random
import itertools

# parse the input 
records = SeqIO.parse(sys.argv[1], "fasta")
sequences = [r.seq for r in records]

k_length = 3


def get_random_str(seq):
    seq_string = str(seq)
    seq_string = seq_string.replace("B", "")
    seq_string = seq_string.replace("X", "")
    seq_string = seq_string.replace("Z", "")
    seq_string = seq_string.replace("U", "")
    seq_string = seq_string.replace("O", "")
    substr_len = random.randrange(20, 41)
    substr = ""
    if len(seq_string) > substr_len + 1:
        index = random.randrange(0, len(seq_string) - substr_len + 1) 
        substr = seq_string[index : (index+substr_len)]
    return substr


def shuffle_natural_peps(natural_peps):
    shuffled = []
    for pep in natural_peps:
        str_var = list(pep)
        random.shuffle(str_var)
        shuffled.append(''.join(str_var))
    return shuffled
    



def produce_peps(seqs):
    natural_peps = []
    random_peps = []

    for seq in seqs:
        random_substr = get_random_str(seq)
        natural_peps.append(random_substr)
        str_var = list(random_substr)
        random.shuffle(str_var)
        random_peps.append(''.join(str_var))
    
    return natural_peps, random_peps

natural_peps, random_peps = produce_peps(sequences)

def create_input_file():
    f = open(sys.argv[2], "w")

    input_length = 20000

    input_natural_peps = []
    input_random_peps = []

    for i in range(input_length):
        random_substr = get_random_str(sequences[i])
        input_natural_peps.append(random_substr)
        
        str_var = list(random_substr)
        random.shuffle(str_var)
        input_random_peps.append(''.join(str_var))
      
    for pep in input_natural_peps:
        f.write(pep)
        f.write("\n")
    
    for pep in input_random_peps:
        f.write(pep)
        f.write("\n")




def create_training_file():
    f = open("train.txt", "w")

    natural_peps = []
    random_peps = []

    for seq in sequences:
        random_substr = get_random_str(seq)
        if random_substr != "":
            natural_peps.append(random_substr)
        
        str_var = list(random_substr)
        random.shuffle(str_var)
        random_peps.append(''.join(str_var))
      
    f.write(str(len(natural_peps)))
    f.write("\n")
    for pep in natural_peps:
        f.write(pep)
        f.write("\n")
    
    for pep in random_peps:
        f.write(pep)
        f.write("\n")



def create_k_mer_freq():
    a_acid_lst = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    comb_lst = []

    for i in range(20):
        for j in range(20):
            for n in range(20):
                lst_var = [a_acid_lst[i], a_acid_lst[j], a_acid_lst[n]]
                comb_lst.append(''.join(lst_var))


    natural_kmer_dict = {}
    random_kmer_dict = {}

    for comb in comb_lst:
        natural_kmer_dict[comb] = 0.0000001
        random_kmer_dict[comb] = 0.0000001


    for pep in natural_peps:
        for i in range(len(pep)-k_length+1):
            natural_kmer_dict[pep[i: i+k_length]] += 1
            
    for pep in random_peps:
        for i in range(len(pep)-k_length+1):
            random_kmer_dict[pep[i: i+k_length]] += 1


    nat_factor = 1.0/ sum(natural_kmer_dict.values())
    for key in natural_kmer_dict.keys():
        natural_kmer_dict[key] = natural_kmer_dict[key]*nat_factor

    ran_factor = 1.0/ sum(random_kmer_dict.values())
    for key in random_kmer_dict.keys():
        random_kmer_dict[key] = random_kmer_dict[key]*ran_factor

    with open("kmer_freq.txt", 'w') as f: 
        for key, value in natural_kmer_dict.items(): 
            f.write('%s:%s\n' % (key, value))

    with open("kmer_freq.txt", 'a') as f: 
        for key, value in random_kmer_dict.items(): 
            f.write('%s:%s\n' % (key, value))

create_input_file()
create_k_mer_freq()
#create_training_file()


