from sklearn import svm
import numpy as np
import tensorflow as tf
from sklearn.metrics import roc_auc_score
import pickle


a_acid_lst = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


natural_pep = []
random_pep = []

natural_vec = []
random_vec = []


trian_file = open('train.txt', 'r')
train_lines = trian_file.readlines()



num_pep = int(train_lines[0].strip())

for i in range(1, num_pep):
    natural_pep.append(train_lines[i].strip()) 

for i in range(num_pep, 2*num_pep+1):
    random_pep.append(train_lines[i].strip()) 



def transform_protein_vector(peps):
  vecs = []
  for pep in peps:
    vec = []
    char_arr = list(pep)
    for char in char_arr:
      vec.append(a_acid_lst.index(char)+1)
    vecs.append(vec)
  
  return vecs


smaple_len = 30000

natural_vecs = transform_protein_vector(natural_pep)[:smaple_len]
random_vecs = transform_protein_vector(random_pep)[:smaple_len]
comb_vecs = natural_vecs + random_vecs

padded_comb_vecs = tf.keras.preprocessing.sequence.pad_sequences(comb_vecs)

train_labels = np.concatenate((np.ones(smaple_len), np.zeros(smaple_len)))

regr = svm.SVR()
regr.fit(padded_comb_vecs, train_labels)


output_filename = 'score2_model.sav'
pickle.dump(regr, open(output_filename, 'wb'))



input_peps = []
input_file = open("input.txt", 'r')
input_lines = input_file.readlines()

   
for line in input_lines:
  input_peps.append(line.strip())


input_vecs = tf.keras.preprocessing.sequence.pad_sequences(transform_protein_vector(input_peps))
input_size = int(len(input_peps) / 2)
input_labels = np.concatenate((np.ones(input_size), np.zeros(input_size)))


predictions = regr.predict(input_vecs)



auroc_score2 = roc_auc_score(input_labels, predictions)
print(auroc_score2)

    



        