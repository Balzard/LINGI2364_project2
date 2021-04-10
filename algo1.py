import copy
# For the Reuters dataset, the negative class will be the file acq.txt and the positive class will be earn.txt. For
# the protein dataset, the negative class will be the file PKA_group15.txt and the positive class will be SRC1521.txt

class Dataset:
    def __init__(self, filepath):
        try:
            data = open(filepath, "r").readlines()
            self.items = list({elem[0]:0 for elem in data if elem[0] != "\n"}.keys())
            self.transactions = []
            i = 1
            while i < len(data):
                if i > 1:
                    i += 1
                if i == len(data):
                    break
                tmp = []
                while data[i] != "\n":
                    tmp.append(data[i][0])
                    i += 1
                tmp = "".join(tmp)
                self.transactions.append(tmp)
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

def spade(filepath, minSupp):
    dataset = Dataset(filepath)
    rep_vert = dataset.to_vertical_representation()
    items = dataset.get_items()
    for item in items:
        supp_tmp = len(set([i[0] for i in rep_vert[item]]))
        if supp_tmp < minSupp:
            continue
        else:
            rep_vert_tmp = copy.deepcopy(rep_vert)
            items_tmp = copy.deepcopy(items)
            items_tmp.remove(item)
            items_tmp.insert(len(items_tmp),item)
            for keys in items_tmp:
                counter1 = 0 # index of first list (item loop)
                counter2 = 0 # index of second list (keys loop)
                # traverse both lists
                while counter1 < len(rep_vert_tmp[item]) and counter2 < len(rep_vert_tmp[keys]):
                    if keys != item:
                        if rep_vert_tmp[keys][counter2][0] == rep_vert_tmp[item][counter1][0]:
                            if rep_vert_tmp[keys][counter2][1] < rep_vert_tmp[item][counter1][1]:
                                rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter2])
                            else:
                                counter2 +=1
                        elif rep_vert_tmp[keys][counter2][0] > rep_vert_tmp[item][counter1][0]:
                            counter1 += 1

                        else:
                            rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter2])
                    else:

                        if rep_vert_tmp[keys][counter1][0] == rep_vert_tmp[keys][counter1 + 1][0] :
                            while rep_vert_tmp[keys][counter1][0] == rep_vert_tmp[keys][counter2][0]:
                                try:
                                    counter2 += 1
                                    rep_vert_tmp[keys][counter2]
                                except:
                                    break
                            rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter1])
                            counter2 -= 1
                            counter1 = counter2
                            if counter1 == len(rep_vert_tmp[keys]) -1:
                                rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter1])

                        else:
                            rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter1])
                            if counter1 == len(rep_vert_tmp[keys]) -1:
                                rep_vert_tmp[keys].remove(rep_vert_tmp[keys][counter1])
                    
        print(item, rep_vert_tmp)

spade("./Datasets/Test/negative.txt",2)

