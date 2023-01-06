import numpy as np
import pandas as pd
import pickle
import os
import Bio
from xgboost import XGBClassifier
from sklearn.multioutput import MultiOutputClassifier
import glob
from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser()                                               

parser.add_argument("--file", "-f", type=str, help='Path to fasta file', required=True)
parser.add_argument("--th1", "-t1", type=float, help='Threshold for Ribosome', default=0.3079)
parser.add_argument("--th2", "-t2", type=float, help='Threshold for Cytosol', default=0.1468)
parser.add_argument("--th3", "-t3", type=float, help='Threshold for Endoplasmic Reticulum', default=0.1156)
parser.add_argument("--th4", "-t4", type=float, help='Threshold for Membrane', default=0.1958)
parser.add_argument("--th5", "-t5", type=float, help='Threshold for Nucleus', default=0.7028)
parser.add_argument("--th6", "-t6", type=float, help='Threshold for Exosome', default=0.9961)
parser.add_argument("--output", "-o", type=str, help='Path to output', required=True)

nf_path = os.path.dirname(__file__)
args = parser.parse_args()

with open(nf_path + '/xgboost_final.pkl', 'rb') as f:
    model = pickle.load(f)

fasta_loc = args.file
th1 = args.th1 
th2 = args.th2 
th3 = args.th3 
th4 = args.th4 
th5 = args.th5 
th6 = args.th6 
output_loc = args.output

exec_string = "python3 " + nf_path + "/Standalone_Nfeature/DNA_Standalone/Nfeature_DNA.py -i " + fasta_loc + " -ft CDK -k 3 -o " + output_loc + "/CDK3.csv"
os.system(exec_string)

print(exec_string)

exec_string = "python3 " + nf_path + "/Standalone_Nfeature/DNA_Standalone/Nfeature_DNA.py -i " + fasta_loc + " -ft RDK -k 4 -o " + output_loc + "/RDK4.csv"
os.system(exec_string)	


# In[103]:


df1 = pd.read_csv(output_loc + "/CDK3.csv")
df1.drop(columns = 'Sequence_ID', inplace=True)
df2 = pd.read_csv(output_loc + "/RDK4.csv")
df2.drop(columns = 'Sequence_ID', inplace=True)
frames = [df1, df2]
df3 = pd.concat(frames, axis=1, ignore_index=True)
df3.shape


y_pred_proba = model.predict_proba(np.array(df3))
y=[]
for i in range(6):
    x = y_pred_proba[i]
    y.append(x[:,1])
y_pred = np.transpose(np.vstack(y))


# In[106]:


motifs_loc = glob.glob(nf_path + "/motifs/*", )
index = [int(j.split('.')[0]) for j in [i.split('_')[-1] for i in motifs_loc]]
motifs_loc_modified = motifs_loc[:]
for i in range(len(index)):
    motifs_loc_modified[index[i]] = motifs_loc[i]

merci_motifs = []

for i in range(len(motifs_loc_modified)):
    with open(motifs_loc_modified[i]) as f:
        merci_motifs.append([line.rstrip() for line in f])    
    
sequences = [str(seq_record.seq) for seq_record in SeqIO.parse(fasta_loc, "fasta")]
identifiers = [seq_record.id for seq_record in SeqIO.parse(fasta_loc, "fasta")]

motif_occurences = []
for i in range(len(merci_motifs)):
    x = merci_motifs[i]
    dummy_list = []
    for j in range(len(x)):
        dummy_list.append([int(merci_motifs[i][j] in k) for k in sequences])

    combined_scores = np.array(dummy_list)
    one_count = np.count_nonzero(combined_scores, axis=0)
    one_count[one_count != 0] = 1
    motif_occurences.append(one_count)


y_pred_motif = y_pred + np.transpose(np.array(motif_occurences))
y_pred_motif[y_pred_motif > 1] = 1 

thresholds = [th1, th2, th3, th4, th5, th6]

prob_df = pd.DataFrame(np.array(y_pred_motif), columns = ['Ribosome', 'Cytosol', 'ER', 'Membrane', 'Nucleus', 'Exosome'])

prob_df_copy = prob_df.copy()
prob_df_copy.insert(loc = 0, column = 'Seq ID', value = identifiers)

for i in range(6):
    x = prob_df.iloc[:,i] 
    prob_df.iloc[:,i][x > thresholds[i]] = 'Yes'
    prob_df.iloc[:,i][x <= thresholds[i]] = 'No'

prob_df.insert(loc = 0, column = 'Seq ID', value = identifiers)
prob_df_copy.to_csv(output_loc +'/final_prob_prediction.csv', header=True, index=False)
prob_df.to_csv(output_loc +'/final_prediction.csv', header=True, index=False)
