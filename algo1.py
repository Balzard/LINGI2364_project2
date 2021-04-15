import copy
import numpy as np
# For the Reuters dataset, the negative class will be the file acq.txt and the positive class will be earn.txt. For
# the protein dataset, the negative class will be the file PKA_group15.txt and the positive class will be SRC1521.txt

class Dataset:
    def __init__(self, filepath):
        try:
            data = open(filepath, "r").readlines()
            self.items = list({elem.split()[0]:0 for elem in data if elem[0] != "\n"}.keys())
            self.transactions = []
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
                except:
                    self.transactions.append(tmp)
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

    def to_vertical_representation(self):
        vet_representation = {i:[] for i in self.items}
        for i in range(len(self.transactions)):
            for j in range(len(self.transactions[i])):
                vet_representation[self.transactions[i][j]].append((i,j))
        return vet_representation


# number transaction starts at 0
# item index in transaction starts at 0


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


def new_explore_branch(rep_vert_pos, rep_vert_neg, branch, counter, k, frequent_seq):
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
        if total_supp_tmp == 0:
            return frequent_seq
        frequent_seq["-".join(branch_tmp)] = [total_supp_tmp, supp_tmp_pos, supp_tmp_neg]
        counter_tmp = counter + 1
        f = new_explore_branch(rep_vert_pos_tmp,rep_vert_neg_tmp, branch_tmp, counter_tmp, k, frequent_seq)
        frequent_seq = merge_two_dicts(frequent_seq,f)

    return frequent_seq


def new_spade(pos_path, neg_path, k):
    assert k > 0
    pos_dataset = Dataset(pos_path)
    neg_dataset = Dataset(neg_path)
    rep_vert_pos = pos_dataset.to_vertical_representation()
    rep_vert_neg = neg_dataset.to_vertical_representation()
    items = sorted(list(set(pos_dataset.get_items()).intersection(neg_dataset.get_items()))) # take items in common
    frequent_seq = {}
    counter = 1
    supports = set()

    for item in items:
        supp_tmp_pos = len(list(dict(rep_vert_pos[item]).items()))
        supp_tmp_neg = len(list(dict(rep_vert_neg[item]).items()))
        total_supp_tmp = supp_tmp_neg + supp_tmp_pos
        frequent_seq[item] = [total_supp_tmp, supp_tmp_pos, supp_tmp_neg] 
    supports = sorted(set([i[0] for i in frequent_seq.values()]))[-k:]

    while len(supports) < k:
        print(items)
        for el in items:
            items_tmp = copy.deepcopy(items)
            items_tmp.remove(el)
            items_tmp.insert(len(items_tmp),el)
            for key in items_tmp:
                branch = [el] + [key]
                rep_vert_pos_tmp = copy.deepcopy(rep_vert_pos)
                rep_vert_neg_tmp = copy.deepcopy(rep_vert_neg)
                rep_vert_pos_tmp = remove_occurences(rep_vert_pos_tmp, el, key)
                rep_vert_neg_tmp = remove_occurences(rep_vert_neg_tmp, el, key) 
                supp_tmp_pos = len(list(dict(rep_vert_pos_tmp[key]).items()))
                supp_tmp_neg = len(list(dict(rep_vert_neg_tmp[key]).items()))
                total_supp_tmp = supp_tmp_neg + supp_tmp_pos
                frequent_seq["-".join(branch)] = [total_supp_tmp, supp_tmp_pos, supp_tmp_neg]
                counter_tmp = counter + 1
                f = new_explore_branch(rep_vert_pos_tmp, rep_vert_neg_tmp, branch, counter_tmp, k, frequent_seq)
                frequent_seq = merge_two_dicts(frequent_seq, f)
                #frequent_seq = dict(list(frequent_seq.items()) + list(new_explore_branch(rep_vert_pos_tmp, rep_vert_neg_tmp, branch, counter_tmp, k, frequent_seq).items())) 
    
        supports = sorted(set([i[0] for i in frequent_seq.values()]))[-k:]
        print(supports)
    
    min_support = supports[0]
    t = []
    #print(frequent_seq)
    for i,j in frequent_seq.items():
        if j[0] >= supports[-k]:
            i = str(i).split("-")
            print(f"{i} {j[1]} {j[2]} {j[0]}")
            t.append(i)
            #t.append(1)
            if j[0] == min_support:
                pass
    return t

    



# t = spade("./Datasets/Test/negative.txt",2)
# print(t)
# print("-------------------------------")
# spade("./Datasets/Test/positive.txt",1)
# spade("./Datasets/Test/negative.txt",1)
#main("./Datasets/Test/positive.txt","./Datasets/Test/negative.txt",2)
#new_spade("./Datasets/Protein/SRC1521.txt", "./Datasets/Protein/PKA_group15.txt",50)
#new_spade("./Datasets/Test/positive.txt","./Datasets/Test/negative.txt",6)
#print(remove_occurences({'B': [(0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 5), (2, 2), (2, 3)], 'A': [(0, 3), (1, 4), (2, 4), (2, 5)],"C": [(0,0),(1,1),(1,3),(2,0),(2,1)]},"A","C"))
# print(remove_occurences({"C":[(0,0),(1,1),(1,3),(2,0),(2,1)]},"C","C"))

# ['B-A-C-A'] 1 3 4

l = [["C"],["C", "A"],["A"],["A", "A"],["B"],["B", "A"],["C", "C"],["C", "C", "A"] ,["A", "C"],["A", "C", "A"],["A", "B"], ["B", "B"],["C", "B"],["C", "B", "A"],["C","B", "B"],["A", "C", "C"] ,["A", "C", "C", "A"] ,["A", "B", "A"],["A", "B", "B"],["C", "C", "B"],
["C", "A", "B"],
["C", "B", "A", "B"],
["A", "A", "A"],
["A", "A", "C"] ,
["A", "A", "C", "A"] ,
["A", "A", "B"],
["A", "B", "C"],
["A", "B", "C", "A"],
["B", "C"] ,
["B", "C", "A"],
["B", "A", "B"] ,
["B", "B", "A"]]
# [C, C, A, A] 
# [C, C, A, B]
# [C, C, B, A] 
# [C, C, B, B] 
# [C, A, A] 
# [C, B, B, A]
# [A, C, C, A, B] 
# [A, C, C, B] 
# [A, C, A, B] 
# [A, C, B]
# [A, C, B, A] 
# [A, C, B, A, B] 
# [A, C, B, B] 
# [A, B, A, B] 
# [B, A, A] 
# [B, B, A, A] ]

r = new_spade("./Datasets/Test/positive.txt","./Datasets/Test/negative.txt",5)

u = [i for i in l if i not in r]
print(u)


