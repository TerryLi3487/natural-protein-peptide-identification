import numpy as np
import sys
import math
from sklearn.metrics import roc_auc_score
import pickle
import tensorflow as tf


k_length = 3 
a_acid_lst = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
input_peps = []

natural_kmer_dict = {}
random_kmer_dict = {}

def transform_protein_vector(peps):
  vecs = []
  for pep in peps:
    vec = []
    char_arr = list(pep)
    for char in char_arr:
      vec.append(a_acid_lst.index(char)+1)
    vecs.append(vec)
  
  return vecs

def calculate_score1(input_peps):
    scores = []
    for pep in input_peps:
        kmer_score_total = 0.0000001
        for i in range(len(pep)-k_length+1):
            kmer = pep[i: i+k_length]
            p = natural_kmer_dict[kmer]
            q = random_kmer_dict[kmer]
            kmer_score_total += math.log2(p / q)  
        scores.append(kmer_score_total)
    return scores

def main():
    kmer_file = open('kmer_freq.txt', 'r')
    kmer_lines = kmer_file.readlines()

    # Strips the newline character
    for i in range(8000):
        key, val = kmer_lines[i].strip().split(":")
        natural_kmer_dict[key] = float(val)

        
    for i in range(8000, 16000):
        key, val = kmer_lines[i].strip().split(":")
        random_kmer_dict[key] = float(val)

    input_file = open(sys.argv[1], 'r')
    input_lines = input_file.readlines()

   
    for line in input_lines:
        input_peps.append(line.strip())


    scores1 = calculate_score1(input_peps)
    model = pickle.load(open("score2_model.sav", 'rb'))
    input_vecs = tf.keras.preprocessing.sequence.pad_sequences(transform_protein_vector(input_peps))
    scores2 = model.predict(input_vecs)

    with open(sys.argv[2], 'a') as f: 
        for i in range(len(input_peps)):
            f.write('%s    %s    %s\n' % (scores1[i], scores2[i], input_peps[i]))


    N = int(len(input_peps) / 2)
    labels = np.concatenate((np.ones(N), np.zeros(N)))
    
    auroc_score1 = roc_auc_score(labels, scores1)
    auroc_score2 = roc_auc_score(labels, scores2)
    print(auroc_score1)
    print(auroc_score2)



if __name__ == "__main__":
    main()
