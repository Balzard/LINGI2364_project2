from datetime import datetime
import numpy as np 

# For the Reuters dataset, the negative class will be the file acq.txt and the positive class will be earn.txt. For
# the protein dataset, the negative class will be the file PKA_group15.txt and the positive class will be SRC1521.txt

dataset = open("./Datasets/Protein/PKA_group15.txt","r")
t = dataset.readlines()