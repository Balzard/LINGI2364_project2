import copy
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

def explore_branch(rep_vert, minSupp, branch, frequent_seq):
    if rep_vert == {}:
        return
    items = list(rep_vert.keys())
    last_item = branch[-1]
    items.remove(last_item)
    items.insert(len(items), last_item)
    rep_vert_tmp = copy.deepcopy(rep_vert)
    for item in items:
        rep_vert_tmp = remove_occurences(rep_vert_tmp, last_item ,item)
        branch_tmp = branch + [item]
        supp_tmp = len(list(dict(rep_vert_tmp[item]).items()))
        if supp_tmp < minSupp:
                del rep_vert_tmp[item]
        else:
            #branch_tmp = branch + [item]
            #print(f"Support of {branch_tmp} is {supp_tmp}")
            frequent_seq["-".join(branch_tmp)] = supp_tmp
            explore_branch(rep_vert_tmp, minSupp, branch_tmp, frequent_seq)
    return frequent_seq


def spade(filepath, minSupp):
    dataset = Dataset(filepath)
    rep_vert = dataset.to_vertical_representation()
    items = dataset.get_items()
    frequent_seq = {}
    for item in items:
        supp_tmp = len(list(dict(rep_vert[item]).items()))
        if supp_tmp < minSupp:
            continue
        else:
            #print(f"Support of {[item]} is {supp_tmp}")
            frequent_seq["".join(item)] = supp_tmp
            rep_vert_tmp = copy.deepcopy(rep_vert)
            items_tmp = copy.deepcopy(items)
            items_tmp.remove(item)
            items_tmp.insert(len(items_tmp),item)
            for keys in items_tmp:
                branch = [item,keys]
                rep_vert_tmp = remove_occurences(rep_vert_tmp, item, keys) 
                if len(list(dict(rep_vert_tmp[keys]).items())) < minSupp:
                    del rep_vert_tmp[keys]
                else:
                    #print(f"Support of {branch} is {len(list(dict(rep_vert_tmp[keys]).items()))}")
                    frequent_seq["-".join(branch)] = supp_tmp
                    d = explore_branch(rep_vert_tmp,minSupp, branch, frequent_seq)
                    frequent_seq = merge_two_dicts(frequent_seq,d)
    return frequent_seq


def main(pos, neg, minSupp):
    positive = spade(pos,1)
    negative = spade(neg, 1)
    c = 0
    for i in positive.keys():
        try:
            t = negative[i] + positive[i]
            tmp = i.split("-")
            if t >= minSupp:
                print(f"{tmp} {positive[i]} {negative[i]} {t}")
                c += 1
            del negative[i]
        except:
            t = positive[i]
            tmp = i.split("-")
            if t >= minSupp:
                print(f"{tmp} {t} {0} {t}")
                c += 1

    if negative != {}:
        for j in negative.keys():
            if negative[j] >= minSupp:
                tmp = j.split("-")
                print(f"{tmp} {0} {negative[j]} {negative[j]}")
                c += 1
    print(c)



# t = spade("./Datasets/Test/negative.txt",2)
# print(t)
# print("-------------------------------")
# spade("./Datasets/Test/positive.txt",2)
main("./Datasets/Test/positive.txt","./Datasets/Test/negative.txt",2)
