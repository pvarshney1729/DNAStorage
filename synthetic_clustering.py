import math
import random
import editdistance
import itertools
import distance
import csv
import sys
import configparser

config_location = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_location)
parameters = config['parameters']
file_locations = config['file_locations']
flags = config['flags']

random.seed()
q_gram_dict = {}
error_rate = parameters['error_rate']
q_gram_call_count = 0
iterable = ['A', 'G', 'C', 'T']
strand_length = int(parameters['strand_length']) # try 110 or 150
threshold = float(parameters['recycle_threshold'])
w = int(parameters['w'])
l = int(parameters['l'])
local_iter = int(parameters['local_iterations'])
global_iter = int(parameters['global_iterations'])
theta_low = int(parameters['theta_low'])
theta_high = int(parameters['theta_high'])
r = int(parameters['r'])
gamma = float(parameters['accuracy_gamma'])
distance_type = "edit"
max_diff = strand_length*(1-threshold)


#Computes the distance (edit/hamming) between 2 strings s1 and s2. Edit distance is O(n^2) while hamming distance computation is O(n) 
def compute_distance(s1, s2, type):
    if type == "edit":
        return editdistance.eval(s1,s2)
    elif type == "hamming":
        distance = sum(c1!= c2 for c1,c2 in zip(s1,s2))
        return distance
    else:
        assert False

#Computes the q_gram representative for a string. Here length = value of q
def q_gram_string(DNA_sequence, length = 3):
    a = itertools.combinations_with_replacement(iterable, length)   #generates all possible combinations from alphabet = iterable of size length. In our case 4^length
    perm = ["".join(ele) for ele in a] #stores them in a list of strings
    #print(perm)
    temp = ''   #equivalent q_gram of DNA_sequence
    for i in perm:  #for each combination check whether it is a substring in the DNA_sequence string. If yes then append 1 else 0.
        if i in DNA_sequence:
            temp += '1'
        else:
            temp += '0'
    return temp #returns the q_gram of the DNA sequence

#generates q=3_gram for blocks of 22 characters in the DNA strand and concatenates them. The resultant is the equivalent q_gram for the DNA strand
def q_gram_sequence(DNA_sequence):
    result = ''
    i = 0
    while i < len(DNA_sequence): #DJORDJE: you should never modify the control variable (i) in the loop like this. You can use while. I'm not sure this is correct. Also, why would you do this in 22-character chunks?
        result += q_gram_string(DNA_sequence[i: i + 20])    #evaluates q=3_gram for consecutive 22 character blocks and appends to results
        i += 20
    return result

def random_string(alphabet = iterable, length = 4):
    return ''.join(random.choices(iterable, k=length))

def initial_clustering(pool, num_strands):
    clusters = [[i] for i in pool]
    return clusters

def in_clusters(clusters, ele): #DJORDJE: what do you use this for?
    for cluster in clusters:
        if ele in cluster: 
            return True
    return False

def clustering(clusters, global_iter, local_iter, w, l):
    new_clusters = clusters
    for i in range(0, global_iter): #in each global iteration
        if(flags.getboolean('print_global_iterations')):
            print("Global iter :", i+1, "#clusters :", len(new_clusters)) #print to screen
        if (i==0 and (int(flags['first_iteration_using_index']==1))):    #for the first iteration create partitions using prefix as hash value
            global_partition = prefix_hashing(new_clusters, w, l)
        else:
            global_partition = global_hashing(new_clusters, w, l)
        global_partition = global_hashing(new_clusters, w, l)
        #print(global_partition)
        new_clusters = []
        global_partition = cluster_partitions(global_partition, w, l)
        for key, val in global_partition.items():
            #print(val)
            for item in val:
                new_clusters.append(item)
        #print(new_clusters)
    return new_clusters

def prefix_hashing(clusters, w, l):
    global_partition = {}   #intialise partitions
    for cluster in clusters:
        #print(cluster)
        representative = random.choice(cluster)
        hash_value = representative[:11]  #DJORDJE: length might need tuning
        if hash_value not in global_partition:
            global_partition[hash_value] = []
        global_partition[hash_value].append(cluster)
    if(flags.getboolean('print_global_partition')):
        count = 1
        for key, cluster in global_partition.items():
            print("partition :", count, "key ", key, "size :", len(cluster))
            count += 1
    #print(global_partition)
    return global_partition
     

def global_hashing(clusters, w, l):
    global_partition = {}
    anchors = []
    num_anchors = int(config['DEFAULT']['global_num_anchors'])
    num_anchor_lists = int(config['DEFAULT']['global_num_anchor_lists'])
    for j in range(0, num_anchor_lists):
        temp = []
        for i in range(0, num_anchors):
            temp.append(random_string(iterable, w))
        anchors.append(temp)
    for cluster in clusters:
        hash_value = ""
        #print(cluster)
        representative = random.choice(cluster)
        #print(representative)
        for j in range(0, num_anchor_lists):
            pos = -1
            for i in range(0, num_anchors):
                if (pos == -1): 
                    pos = representative.find(anchors[j][i], 0, len(representative) - (w+l)) #DJORDJE: should be double checked
                else:
                    break
            m_prime = pos+w+l
            hash_value += representative[pos:m_prime]
        #print(anchor, hash_value)
        if hash_value not in global_partition:
            global_partition[hash_value] = []
        global_partition[hash_value].append(cluster)
    if(flags.getboolean('print_global_partition')):
        count = 1
        for key, cluster in global_partition.items():
            print("partition :", count, "key ", key, "size :", len(cluster))
            count += 1
    #print(global_partition)
    return global_partition


def local_clustering(clusters, w, l):
    global q_gram_call_count
    new_clusters = clusters
    freq_representative = int(config['DEFAULT']['frequency_of_choosing_new_random_representative'])
    sorting_or_pairwise = flags.getboolean('sorting_or_pairwise')

    if(sorting_or_pairwise):
        if(len(new_clusters) >= 2):
            for k in range(0, local_iter):
                if(k%freq_representative == 0):
                    anchors1 = []
                    anchors2 = []
                    for i in range(0, 20):
                        anchors2.append(random_string(iterable, w))
                        anchors1.append(random_string(iterable, w))
                    hash_values = []
                    for cluster in new_clusters:
                        hash_value = ""
                        #print(cluster)
                        representative = random.choice(cluster)
                        #print(representative)
                        pos1 = -1
                        pos2 = -1
                        for i in range(0, 20):
                            if (pos1 == -1 and pos2 == -1):
                                pos1 = representative.find(anchors1[i], 0, len(representative) - (w+l)) #DJORDJE: you need to tripple check the semantics of find regarding the end of your search. For example, I don't think this would find an anchor that sits between positions l+w and l (from the end)
                                pos2 = representative.find(anchors2[i], 0, len(representative) - (w+l))
                            elif (pos1 == -1):
                                pos1 = representative.find(anchors1[i], 0, len(representative) - (w+l))
                            elif (pos2 == -1):
                                pos2 = representative.find(anchors2[i], 0, len(representative) - (w+l))
                            else:
                                break
                #DJORDJE: Did you test the code above?
                        m_prime1 = pos1+w+l
                        m_prime2 = pos2+w+l

                #DJORDJE: There is a chance that a) no anchor is found in one or two of the lists and b) the anchors found in the two lists can be identical. I don't think these cases are handled properly 
                        hash_value = representative[pos1:m_prime1] + representative[pos2:m_prime2]
                        #print(anchor, hash_value)
                        hash_values.append([hash_value, representative, cluster])

                    hash_values.sort(key= lambda x:x[0]) 
                    #print(k, hash_values)
                    i = 0
                    n = len(hash_values)                                                                 
                    while i < n-1:
                        #print(i, n)
                        representative1 = hash_values[i][1]
                        j = i + 1
                        representative2 = hash_values[j][1]
                        q_gram_call_count += 2
                        hamming_distance = compute_distance(q_gram_dict[representative1], q_gram_dict[representative2], "hamming")
                        if(hamming_distance <= theta_low):
                            hash_values[i][2] += hash_values[j][2]
                            hash_values[i][1] = random.choice(hash_values[i][2])  #DJORDJE: Perhaps we should pick a new representative after merging
                            del hash_values[j]
                            n = n-1
                            #print(n)
                        elif(hamming_distance <= theta_high):
                            edit_distance = compute_distance(representative1, representative2, "edit")
                            if (edit_distance <= r):
                                hash_values[i][2] += hash_values[j][2] 
                                hash_values[i][1] = random.choice(hash_values[i][2])
                                del hash_values[j]
                                n = n-1
                                #print(n)
                            else:
                                i+=1
                        else:
                            i += 1
                    new_clusters = []
                    i = 0
                    while i < len(hash_values):
                        new_clusters.append(hash_values[i][2])
                        i += 1
                
                else:
                    anchors1 = []
                    anchors2 = []
                    for i in range(0, 20):
                        anchors2.append(random_string(iterable, w))
                        anchors1.append(random_string(iterable, w))
                    j = 0
                    n = len(hash_values)
                    while j < n:
                        #print(hash_values[i])
                        hash_value = ""
                        #print(cluster)
                        representative = random.choice(hash_values[j][1])
                        #print(representative)
                        pos1 = -1
                        pos2 = -1
                        for i in range(0, 20):
                            if (pos1 == -1 and pos2 == -1):
                                pos1 = representative.find(anchors1[i], 0, len(representative) - (w+l)) #DJORDJE: you need to tripple check the semantics of find regarding the end of your search. For example, I don't think this would find an anchor that sits between positions l+w and l (from the end)
                                pos2 = representative.find(anchors2[i], 0, len(representative) - (w+l))
                            elif (pos1 == -1):
                                pos1 = representative.find(anchors1[i], 0, len(representative) - (w+l))
                            elif (pos2 == -1):
                                pos2 = representative.find(anchors2[i], 0, len(representative) - (w+l))
                            else:
                                break
                        m_prime1 = pos1+w+l
                        m_prime2 = pos2+w+l

                        hash_value = representative[pos1:m_prime1] + representative[pos2:m_prime2]
                        hash_values[j][0] = hash_value
                        j += 1

                    hash_values.sort(key= lambda x:x[0])
                    i = 0
                    n = len(hash_values) #DJORDJE: Shouldn't you look for new hashes?
                    while i < n-1:
                        representative1 = hash_values[i][1]
                        j = i + 1
                        representative2 = hash_values[j][1]
                        q_gram_call_count += 2
                        hamming_distance = compute_distance(q_gram_dict[representative1], q_gram_dict[representative2], "hamming")
                        if(hamming_distance <= theta_low):
                            hash_values[i][2] += hash_values[j][2]
                            hash_values[i][1] = random.choice(hash_values[i][2])
                            del hash_values[j]
                            n = n-1
                            #print(n)
                        elif(hamming_distance <= theta_high):
                            edit_distance = compute_distance(representative1, representative2, "edit")
                            if (edit_distance <= r):
                                hash_values[i][2] += hash_values[j][2]
                                hash_values[i][1] = random.choice(hash_values[i][2])
                                del hash_values[j]
                                n = n-1
                                #print(n)
                            else:
                                i+=1
                        else:
                            i += 1
                    new_clusters = []
                    i = 0
                    while i < len(hash_values):
                        new_clusters.append(hash_values[i][2])
                        i += 1
    else:
        if(len(new_clusters) >= 2):
            for k in range(0, local_iter):
                if(k%freq_representative == 0):
                    num_anchors = int(config['DEFAULT']['local_num_anchors'])
                    num_anchor_lists = int(config['DEFAULT']['local_num_anchor_lists'])
                    anchors = []
                    for j in range(0, num_anchor_lists):
                        temp = []
                        for i in range(0, num_anchors):
                            temp.append(random_string(iterable, w))
                        anchors.append(temp)
                    hash_values = []
                    unique_hash_values = []
                    num = 0
                    for cluster in new_clusters:
                        hash_value = ""
                        #print(cluster)
                        representative = random.choice(cluster)
                        #print(representative)
                        pos1 = -1
                        pos2 = -1
                        for j in range(0, num_anchor_lists):
                            pos = -1
                            for i in range(0, num_anchors):
                                if (pos == -1): 
                                    pos = representative.find(anchors[j][i], 0, len(representative) - (w+l)) #DJORDJE: should be double checked
                                else:
                                    break
                            m_prime = pos+w+l
                            hash_value += representative[pos:m_prime]
                        if hash_value not in unique_hash_values:
                            unique_hash_values.append(hash_value)
                        #print(anchor, hash_value)
                        hash_values.append([hash_value, representative])
                        num +=1

                    unique_hash_values.sort()
                    i = 0
                    n = len(new_clusters)
                    j = 0
                    while i < n:
                        hash_value1 = hash_values[i][0]
                        representative1 = hash_values[i][1]
                        j = i + 1
                        while j < n:
                            hash_value2 = hash_values[j][0]
                            representative2 = hash_values[j][1]
                            if(hash_value1==hash_value2):
                                q_gram_call_count += 2
                                hamming_distance = compute_distance(q_gram_dict[representative1], q_gram_dict[representative2], "hamming")
                                if(hamming_distance <= theta_low):
                                    new_clusters[i] += new_clusters[j]
                                    hash_values[i][1] = random.choice(new_clusters[i])
                                    del new_clusters[j]
                                    del hash_values[j]
                                    n = n-1
                                    #print(n)
                                else:
                                    edit_distance = compute_distance(representative1, representative2, "edit")
                                    if (hamming_distance <= theta_high and edit_distance <= r):
                                        new_clusters[i] += new_clusters[j]
                                        hash_values[i][1] = random.choice(new_clusters[i])
                                        del new_clusters[j]
                                        del hash_values[j]
                                        n = n-1
                                        #print(n)
                                    else:
                                        j += 1
                            else:
                                j = j+1
                        i = i+1
                else:
                    num_anchors = int(config['DEFAULT']['local_num_anchors'])
                    num_anchor_lists = int(config['DEFAULT']['local_num_anchor_lists'])
                    anchors = []
                    for j in range(0, num_anchor_lists):
                        temp = []
                        for i in range(0, num_anchors):
                            temp.append(random_string(iterable, w))
                        anchors.append(temp)
                    num = 0
                    for m in range(len(hash_values)):
                        hash_value = ""
                        #print(cluster)
                        representative = hash_values[m][1]
                        #print(representative)
                        pos1 = -1
                        pos2 = -1
                        for j in range(0, num_anchor_lists):
                            pos = -1
                            for i in range(0, num_anchors):
                                if (pos == -1): 
                                    pos = representative.find(anchors[j][i], 0, len(representative) - (w+l)) #DJORDJE: should be double checked
                                else:
                                    break
                            m_prime = pos+w+l
                            hash_value += representative[pos:m_prime]
                        #print(anchor, hash_value)
                        hash_values[m][0] = hash_value
                        num +=1
                    i = 0
                    n = len(new_clusters)
                    j = 0
                    while i < n:
                        hash_value1 = hash_values[i][0]
                        representative1 = hash_values[i][1]
                        j = i + 1
                        while j < n:
                            hash_value2 = hash_values[j][0]
                            representative2 = hash_values[j][1]
                            if(hash_value1==hash_value2):
                                q_gram_call_count += 2
                                hamming_distance = compute_distance(q_gram_dict[representative1], q_gram_dict[representative2], "hamming")
                                if(hamming_distance <= theta_low):
                                    new_clusters[i] += new_clusters[j]
                                    hash_values[i][1] = random.choice(new_clusters[i])
                                    del new_clusters[j]
                                    del hash_values[j]
                                    n = n-1
                                    #print(n)
                                else:
                                    edit_distance = compute_distance(representative1, representative2, "edit")
                                    if (hamming_distance <= theta_high and edit_distance <= r):
                                        new_clusters[i] += new_clusters[j]
                                        hash_values[i][1] = random.choice(new_clusters[i])
                                        del new_clusters[j]
                                        del hash_values[j]
                                        n = n-1
                                        #print(n)
                                    else:
                                        j += 1
                            else:
                                j = j+1
                        i = i+1
    return new_clusters



def cluster_partitions(global_partition, w, l):
    count = 1
    for key, clusters in global_partition.items():
        global_partition[key] = local_clustering(clusters, w, l)
        #print(count)
        count += 1
    return global_partition

def underlying_clustering(pool, actual_num_strands, multiplier):
    clusters = []   #initialise the underlying clusters
    temp = []   #intialise an empty cluster
    for i in range(0, len(pool)):   #read individual strands from the pool
        temp.append(pool[i])    #append the strand to the current cluster
        if((i+1)%10 == 0):  #if we have reached a multiple of 10 (i.e. block of 10 strands is over)
            j = 0
            clusters.append(temp)   #append the cluster to the underlying_clusters list
            temp = []   #start afresh with an empty cluster
    return clusters

#for each strand 'ele' in pool, stores its q_gram representative in q_gram_dict[ele]. This leads to creation of a dictionary with (strand, equivalent q_gram) pairs
def create_q_gram_dict(pool):
    for ele in pool:
        q_gram_dict[ele] = q_gram_sequence(ele)


def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3

def accuracy(clusters, underlying_clusters):
    count = 0
    accuracy = 0
    for i in range(0, len(underlying_clusters)):
        i_count = 0
        for j in range(0, len(clusters)):
            temp = intersection(clusters[j], underlying_clusters[i])
            if ((temp == clusters[j]) and (len(temp)/len(underlying_clusters[i])) >= gamma) :
                i_count +=1
                break
                #print(i, j)
        #print(i, i_count)
        count += i_count
    accuracy = count/len(underlying_clusters)
    
    return accuracy


def main():
    global gamma
    actual_num_strands = parameters['underlying_strands']   #no. of underlying strands

    #Read the pool of strands from txt file (which has each line representing one strand)
    pool = []   
    with open(file_locations['input_pool'], "r") as f:
        pool = f.read().splitlines()

    #create the q_gram dictionary for each strand in the pool
    create_q_gram_dict(pool)
    
    #stores the underlying clusters from the pool. The txt file generated stores the strands belonging to the same cluster together.
    #i.e. each block of 10 consecutive strands is a cluster.
    underlying_clusters = underlying_clustering(pool, actual_num_strands, 10)

    #print(underlying_clusters[0], underlying_clusters[1],len(underlying_clusters))

    clusters_paper = initial_clustering(pool, actual_num_strands)
    #print(clusters)

    #Implement the algorith with underlying clusters in hand
    clusters = clustering(clusters_paper, global_iter, local_iter, w ,l)
    
    print("Error rate:", error_rate)
    print("Paper clusters:", len(clusters))
    print("Underlying Clusters: ", len(underlying_clusters))
    with open(file_locations['underlying_clusters'], "w") as g:
        wr = csv.writer(g)
        wr.writerows(underlying_clusters)
    
    with open(file_locations["output_clusters"], "w") as g:
        wr = csv.writer(g)
        wr.writerows(clusters)
    
    print("Accuracy (gamma=", gamma,") :", accuracy(clusters, underlying_clusters))
    
if __name__=="__main__":
    main()
