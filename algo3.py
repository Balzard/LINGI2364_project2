#!/usr/bin/env python3
import sys
import copy
import numpy as np
import math
# For the Reuters dataset, the negative class will be the file acq.txt and the positive class will be earn.txt. For
# the protein dataset, the negative class will be the file PKA_group15.txt and the positive class will be SRC1521.txt

class Dataset:
    def __init__(self, filepath):
        try:
            data = open(filepath, "r").readlines()
            self.items = list({elem.split()[0]:0 for elem in data if elem[0] != "\n"}.keys())
            self.transactions = []
            self.longest_trans = 0
            i = 0
            while i < len(data):
                if data[i] == "\n":
                    i += 1
                # if i == len(data) :
                #     break
                tmp = []
                try:
                    while data[i] != "\n":
                        tmp.append(data[i].split()[0])
                        i += 1
                    self.transactions.append(tmp)
                    if len(tmp) > self.longest_trans:
                        self.longest_trans = len(tmp)
                except:
                    self.transactions.append(tmp)
                    if len(tmp) > self.longest_trans:
                        self.longest_trans = len(tmp)
                    break
            
        except IOError as e :
            print("Unable to read dataset file\n" + e)

    def get_items(self):
        return sorted(self.items)

    def get_transactions(self):
        return self.transactions

    def get_transaction(self, i):
        return self.transactions[i]

    def get_nb_transactions(self):
        return len(self.transactions)

    def get_longest_transaction(self):
        return self.longest_trans

    def to_vertical_representation(self):
        vet_representation = {i:[] for i in self.items}
        for i in range(len(self.transactions)):
            for j in range(len(self.transactions[i])):
                vet_representation[self.transactions[i][j]].append((i,j))
        return vet_representation


# number transaction starts at 0
# item index in transaction starts at 0


def wracc(p, n, px, nx):
    return round(((p/(p+n))*(n/(p+n)))*((px/p)-(nx/n)),5)


def remove_occurences(rep_vert, index1, index2):
    """
    Index2 and index1 are items
    Remove element from index2 list that cannot be mapped to index1 list.

    Parameters:
        rep_vert : dictionary where each key is an item and the value is its vertical representation
        index1 (int): dictionary key, item of dataset
        index2 (int): dictionary key, item of dataset

    Returns: 
        rep_vert edited
    """

    if rep_vert[index1] == []:
        rep_vert[index2] = []
        return rep_vert
    
    counter1 = 0 # index of first list (item loop)
    counter2 = 0 # index of second list (keys loop)
    # traverse both lists
    while counter1 < len(rep_vert[index1]) and counter2 < len(rep_vert[index2]):
        if index2 != index1:
            if rep_vert[index2][counter2][0] == rep_vert[index1][counter1][0]:
                if rep_vert[index2][counter2][1] < rep_vert[index1][counter1][1]:
                    rep_vert[index2].remove(rep_vert[index2][counter2])
                else:
                    counter2 +=1
            elif rep_vert[index2][counter2][0] > rep_vert[index1][counter1][0]:
                try:
                    counter1 += 1
                    rep_vert[index1][counter1]
                except:
                    del rep_vert[index2][-(len(rep_vert[index2])-counter2):]

            else:
                rep_vert[index2].remove(rep_vert[index2][counter2])
        else:

            try:
                if rep_vert[index2][counter1][0] == rep_vert[index2][counter1 + 1][0] :
                    while rep_vert[index2][counter1][0] == rep_vert[index2][counter2][0]:
                        try:
                            counter2 += 1
                            rep_vert[index2][counter2]
                        except:
                            break
                    rep_vert[index2].remove(rep_vert[index2][counter1])
                    counter2 -= 1
                    counter1 = counter2
                    if counter1 == len(rep_vert[index2]) -1:
                        rep_vert[index2].remove(rep_vert[index2][counter1])
                
                else:
                    rep_vert[index2].remove(rep_vert[index2][counter1])
                    if counter1 == len(rep_vert[index2]) -1:
                        rep_vert[index2].remove(rep_vert[index2][counter1])
            
            except IndexError:
                rep_vert[index2].remove(rep_vert[index2][counter1])
                    
    return rep_vert

# from stackoverflow
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


def explore_branch(rep_vert_pos, rep_vert_neg, branch, counter, k, frequent_seq, nb_trans_pos, nb_trans_neg):
    """
    Explore in depth a branch sequence

    Parameters:
        rep_vert_pos (dict): vertical representation of positive dataset
        rep_vert_neg (dict): vertical representation of negative dataset
        branch (list): current sequence visited in the search
        counter (int): level to stop in the depth first search
        k (int): number of k most frequent patterns
        frequent_seq (dict): dict where each key is a sequence and its value is an array [total_supp, supp_pos, supp__neg]

    Returns:
        frequent_seq (dict) edited
    """
    if counter == k:
        return frequent_seq
    items = sorted(list(set(rep_vert_pos.keys()).union(list(rep_vert_neg.keys()))))
    last_item = branch[-1]
    items.remove(last_item)
    items.insert(len(items), last_item)
    for item in items:
        rep_vert_pos_tmp = copy.deepcopy(rep_vert_pos)
        rep_vert_neg_tmp = copy.deepcopy(rep_vert_neg)
        rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, last_item, item)
        rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, last_item, item)
        branch_tmp = branch + [item]
        supp_tmp_pos = len(list(dict(rep_vert_pos_tmp[item]).items()))
        supp_tmp_neg = len(list(dict(rep_vert_neg_tmp[item]).items()))
        total_supp_tmp = supp_tmp_neg + supp_tmp_pos
        if total_supp_tmp > 0 and supp_tmp_pos > 0:
            w = wracc(nb_trans_pos, nb_trans_neg, supp_tmp_pos, supp_tmp_neg)
            frequent_seq["-".join(branch_tmp)] = [w, supp_tmp_pos, supp_tmp_neg]
            counter_tmp = counter + 1
            f = explore_branch(rep_vert_pos_tmp, rep_vert_neg_tmp, branch_tmp, counter_tmp, k, frequent_seq,nb_trans_pos, nb_trans_neg)
            frequent_seq = merge_two_dicts(frequent_seq, f)
        try:
            if frequent_seq[branch][1] == supp_tmp_pos and frequent_seq[branch][2] == supp_tmp_neg:
                del frequent_seq[branch]
        except:
            pass
                   

    return frequent_seq


def spade(pos_path, neg_path, k):
    """
    Returns the k most frequent patterns in both datasets

    Parameters:
        pos_path (string): path to positive dataset
        neg_path (string): path to negative dataset
        k (int): number of most frequent patterns to find

    Returns:
        None
    """
    assert k > 0
    pos_dataset = Dataset(pos_path)
    neg_dataset = Dataset(neg_path)
    rep_vert_pos = pos_dataset.to_vertical_representation()
    rep_vert_neg = neg_dataset.to_vertical_representation()
    items = sorted(list(set(pos_dataset.get_items()).intersection(neg_dataset.get_items()))) # take items in common
    items_tmp = copy.deepcopy(items)
    frequent_seq = {}
    counter = 0
    supports = set()
    k_tmp = 2
    nb_trans_pos = pos_dataset.get_nb_transactions() - 1
    nb_trans_neg = neg_dataset.get_nb_transactions() - 1
    
    for item in items:
        supp_tmp_pos = len(list(dict(rep_vert_pos[item]).items()))
        supp_tmp_neg = len(list(dict(rep_vert_neg[item]).items()))
        total_supp_tmp = supp_tmp_neg + supp_tmp_pos
        if total_supp_tmp > 0 and supp_tmp_pos > 0:
            w = wracc(nb_trans_pos, nb_trans_neg, supp_tmp_pos, supp_tmp_neg)
            frequent_seq[item] = [w, supp_tmp_pos, supp_tmp_neg] 
    supports = sorted(set([i[0] for i in frequent_seq.values()]))[-k:]

    while True:
        for el in items:
            if isinstance(el, str):
                items_tmp.remove(el)
                items_tmp.insert(len(items_tmp),el)
            else:
                items_tmp.remove(el[-1])
                items_tmp.insert(len(items_tmp),el[-1])
            for key in items_tmp:
                if len(el) == 1 or isinstance(el, str):
                    branch = [el] + [key] 
                else:
                    branch = el + [key]

                if len(branch) == 2:
                    rep_vert_pos_tmp = copy.deepcopy(rep_vert_pos)
                    rep_vert_neg_tmp = copy.deepcopy(rep_vert_neg)
                    counter_tmp = 2
                else:
                    rep_vert_pos_tmp = copy.deepcopy(rep_vert_pos)
                    rep_vert_neg_tmp = copy.deepcopy(rep_vert_neg)
                    counter_tmp = 1
                    for i in range(len(el)-1):
                        if isinstance(el,str):
                            pass
                        else:
                            rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, branch[i], branch[i+1])
                            rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, branch[i], branch[i+1]) 
                if isinstance(el,str):
                    rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, el, key)
                    rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, el, key)
                else:
                    rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, el[-1], key)
                    rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, el[-1], key) 
                supp_tmp_pos = len(list(dict(rep_vert_pos_tmp[key]).items()))
                supp_tmp_neg = len(list(dict(rep_vert_neg_tmp[key]).items()))
                total_supp_tmp = supp_tmp_neg + supp_tmp_pos
                if total_supp_tmp > 0 and supp_tmp_pos > 0:
                    w = wracc(nb_trans_pos, nb_trans_neg, supp_tmp_pos, supp_tmp_neg)
                    frequent_seq["-".join(branch)] = [w, supp_tmp_pos, supp_tmp_neg]
                    f = explore_branch(rep_vert_pos_tmp, rep_vert_neg_tmp, branch, counter_tmp, k_tmp, frequent_seq, nb_trans_pos, nb_trans_neg)
                    frequent_seq = merge_two_dicts(frequent_seq, f)

                try:
                    if frequent_seq[el][1] == supp_tmp_pos and frequent_seq[el][2] == supp_tmp_neg:
                        print("here"+str(el))
                        del frequent_seq[el]
                except:
                    pass
                 

        if counter == 0:
            counter = 2
        else:
            counter = k_tmp + counter

        supports = sorted(set([i[0] for i in frequent_seq.values()]))[-k:]

        if len(supports) < k:
            items = []
            for i in frequent_seq.keys():
                i_tmp = str(i).split("-")
                if len(i_tmp) == k_tmp:
                    items.append(i_tmp)

        if len(supports) >= k:
            break
    
    min_support = supports[0]
    items = []
    for i,j in list(frequent_seq.items()):
        r_tmp = str(i).split("-")
        if len(r_tmp) == counter:
            items.append(r_tmp)

    while items != []:
        l_tmp = []
        for j_tmp in items:
                items_tmp.remove(j_tmp[-1])
                items_tmp.insert(len(items_tmp),j_tmp[-1])
                for key in items_tmp:
                    rep_vert_pos_tmp = copy.deepcopy(rep_vert_pos)
                    rep_vert_neg_tmp = copy.deepcopy(rep_vert_neg)
                    for i in range(len(j_tmp)-1):
                        rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, j_tmp[i], j_tmp[i+1])
                        rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, j_tmp[i], j_tmp[i+1]) 
                    rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, j_tmp[-1], key)
                    rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, j_tmp[-1], key) 
                    supp_tmp_pos = len(list(dict(rep_vert_pos_tmp[key]).items()))
                    supp_tmp_neg = len(list(dict(rep_vert_neg_tmp[key]).items()))
                    total_supp_tmp = supp_tmp_neg + supp_tmp_pos
                    tmp = j_tmp + [key]
                    w = wracc(nb_trans_pos, nb_trans_neg, supp_tmp_pos, supp_tmp_neg)
                    if total_supp_tmp > 0:
                        items.append(tmp)
                        frequent_seq["-".join(tmp)] = [w, supp_tmp_pos, supp_tmp_neg]
                        if w > min_support:
                            supports.pop(0)
                            supports.append(w)
                            supports = sorted(set([i[0] for i in frequent_seq.values()]))[-k:]
                            min_support = supports[0]
                    try:
                        if frequent_seq[j_tmp][1] == supp_tmp_pos and frequent_seq[j_tmp][2] == supp_tmp_neg:
                            del frequent_seq[j_tmp]
                    except:
                        pass
                        
    
                l_tmp.append(j_tmp)
        counter += 1
        items = [x for x in items if x not in l_tmp]
      
    t = 0
    for i,j in frequent_seq.items():
        i2 = str(i).split("-")
        if j[0] >= min_support:
            print(f"[{', '.join(map(str, i2))}] {j[1]} {j[2]} {j[0]}")
            t+=1
    print(t)
        
    

def main():
    pos_filepath = sys.argv[1] # filepath to positive class file
    neg_filepath = sys.argv[2] # filepath to negative class file
    k = int(sys.argv[3])
    spade(pos_filepath, neg_filepath, k)


if __name__ == "__main__":
    #spade("./Datasets/Protein/SRC1521.txt","./Datasets/Protein/PKA_group15.txt",10)
    spade("./Datasets/Test/positive.txt","./Datasets/Test/negative.txt",6)
    #main()




