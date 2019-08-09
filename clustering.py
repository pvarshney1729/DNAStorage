import math
import random
from bitstring import BitArray
from seq_codec import encodePattern, decodePattern
import editdistance
import itertools
import distance
import pickle
import csv
from PIL import Image
import multiprocessing as mp
import configparser
import sys
import os
import time
start_time = time.time()

config_location = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_location)
parameters = config['parameters']
file_locations = config['file_locations']
flags = config['flags']

actual_num_strands = 0
manager = mp.Manager()
random.seed()
xor_string = []
q_gram_dict = {}
q_gram_call_count = 0
iterable = ['A', 'G', 'C', 'T']
strand_length = int(parameters['strand_length']) # try 110 or 150
error_rate = float(parameters['error_rate'])
threshold = float(parameters['recycle_threshold'])
w = int(parameters['w'])
l = int(parameters['l'])
gamma = float(parameters['accuracy_gamma'])
local_iter = int(parameters['local_iterations'])
global_iter = int(parameters['global_iterations'])
theta_low = int(parameters['theta_low'])
theta_high = int(parameters['theta_high'])
r = int(parameters['r'])
distance_type = "edit"
file_path = file_locations['input_img']
deviation_fraction = float(parameters['deviation_from_cluster'])
wrong_bases_limit = int(parameters['wrong_bases_limit'])
min_cluster_size = int(parameters['min_cluster_size'])
max_edit_limit = int(parameters['max_edit_dist'])
max_diff = strand_length*(1-threshold)

def decode(DNA_array, file_content, k, I, bin_len):
    # decode nucleotides to ternary
    decoded_data = [None] * k  # {decimal_index: binary_data}

    for i in range(0, k):
        ele = DNA_array[i]
        decoded_ternary_idx = ''
        decoded_ternary = ''
        codePost = 'A'
        for j in range(0, I):
            decoded_bit = decodePattern[(codePost, ele[j])]
            decoded_ternary_idx += decoded_bit
            codePost = ele[j]

        codePost = 'A'
        for j in range(I, strand_length):
            decoded_bit = decodePattern[(codePost, ele[j])]
            decoded_ternary += decoded_bit
            codePost = ele[j]

        decimal_index = int(decoded_ternary_idx, 3)
        binary_data = knary_fix_length(int(decoded_ternary, 3), 2, bin_len)
        decoded_data[decimal_index] = binary_data

    # recover the last binary data
    strip_len = k * bin_len - len(file_content)
    decoded_data[-1] = decoded_data[-1][strip_len:]

    bin_data = ''
    for item in decoded_data:
        bin_data += item

    return bin_data

def decodeStrand(DNA_array):
    # decode nucleotides to ternary
    decoded_ternary_idx = ""
    codePost = 'A'
    for j in range(0, len(DNA_array)):
        decoded_bit = decodePattern[(codePost, DNA_array[j])]
        decoded_ternary_idx += decoded_bit
        codePost = DNA_array[j]

    return decoded_ternary_idx

def sxor(s1,s2):    
    return ''.join(chr(int(a) ^ int(b) + 48 ) for a,b in zip(s1,s2))

def read_binary(file_path):
    global xor_string
    f = open(file_path, mode='rb')
    file_content = BitArray(f.read())
    file_content = file_content.bin
    #print(file_content)
    xor_string = ''.join(random.choices(['0','1'], k=len(file_content)))
    #print(xor_string)
    file_content = sxor(file_content, xor_string)
    #print(type(file_content))
    print('Binary length of file = ', len(file_content))
    return file_content

def strands_compute(strand_length, file_length):
     num_strands = 0
     data_length = strand_length    #in ternary
     index_length = 1             # in ternary

     while True:
         data_length_binary = int(data_length*math.log(3,2))
         num_strands = math.ceil(file_length/data_length_binary)
         index_length = math.ceil(math.log(num_strands,3))
         if index_length + data_length == strand_length:
            break
         data_length = data_length - 1
    
     return num_strands, data_length, index_length, data_length_binary

def knary_fix_length(n, base, size):
    if n==0:
        return '0'*size
    ans = []
    while n:
        n, r = divmod(n, base)
        ans.append(str(r))
    if len(ans)<size :
        ans += '0'*(size-len(ans))
    return ''.join(reversed(ans))

def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i + n]

def encode(file_content, num_strands, data_length, index_length, data_length_binary):
    #binary_sequence = list(chunks(file_content, data_length_binary))
    binary_sequence = [file_content[i:i+data_length_binary] for i in range(0,len(file_content),data_length_binary)]
    #print(binary_sequence)
    DNA_sequence = []
    for i in range(0, len(binary_sequence)):
        ele = binary_sequence[i]
        prev = 'A'
        ternary_data = knary_fix_length(int(ele,2), 3, data_length)
        ternary_index = knary_fix_length(i, 3, index_length)

        temp_DNA = ''
        for j in ternary_index:
            prev = encodePattern[(prev, int(j))]
            temp_DNA += prev
        prev = 'A'
        for j in ternary_data:
            prev = encodePattern[(prev, int(j))]
            temp_DNA += prev

        #print(temp_DNA)

        DNA_sequence.append(temp_DNA)
    print(num_strands, index_length, data_length, num_strands*strand_length)
    return DNA_sequence

bases = ['A','C','G','T']
choices = [('sub',100*error_rate),('in',100*error_rate),('del',100*error_rate),('skip',300-300*error_rate)]
total = sum(w for c, w in choices)

def weighted_choice(choices):
   r = random.uniform(0, total)
   upto = 0
   for c, w in choices:
      if upto + w >= r:
         return c
      upto += w
   assert False, "Shouldn't get here"

def add_noise(strand):
    out = []
    for c in strand:
        error = weighted_choice(choices)
        if 'sub' == error:
            out.append(random.choice(bases))
        elif 'in' == error:
            out.append(random.choice(bases))
            out.append(c)
        elif 'skip' == error:
            out.append(c)
        # else delete
    return ''.join(out)

def generate_strand_error(DNA_sequence):
    multiplier = int(parameters['pool_multiplier'])
    pool = DNA_sequence.copy()
    pool *= multiplier
    
    if(flags.getboolean('our_errorgen_or_cyrus')):
        pool_size = len(pool) * strand_length
        insert_pos = random.sample(range(pool_size), int(error_rate * pool_size))
        insert_pos.sort(reverse=True)

        # 0: insert  1: delete
        error_type = [random.randint(0, 2) for i in range(int(error_rate * pool_size))]

        insert_val = [random.randint(0, 2) for i in range(int(error_rate * pool_size))]

        error_sequence = zip(insert_pos, error_type, insert_val)

        for i,  comb in enumerate(error_sequence):

            pos = comb[0]
            error = comb[1]
            val = comb[2]
            strand_idx = int(pos / strand_length)
            base_idx = pos % strand_length

            current_base = pool[strand_idx][base_idx]
            current_strand = list(pool[strand_idx])

            if error == 0:
                current_strand[base_idx] = encodePattern[(current_base, val)]
            elif error ==1:
                del current_strand[base_idx]
            elif error == 2:
                current_strand[base_idx] = random.choice(iterable)
            else:
                assert False
            #print (len(current_strand))
            pool[strand_idx] = "".join(current_strand)

    else:
        for i in range(len(pool)):
                #print(replicate, compute_distance(strand, alpha, "edit"))
            pool[i] = add_noise(pool[i])
    return pool

def decode_strand(DNA_sequence):
    decoded_ternary = ""
    prev = 'A'
    for  i in range(len(DNA_sequence)):
        decoded_base = decodePattern[(prev, DNA_sequence[i])]
        decoded_ternary += decoded_base
        prev = DNA_sequence[i]

    return decoded_ternary

def compute_distance(s1, s2, type):
    if type == "edit":
        return editdistance.eval(s1,s2)
    elif type == "hamming":
        distance = sum(c1!= c2 for c1,c2 in zip(s1,s2))
        return distance
    else:
        assert False

a = itertools.combinations_with_replacement(iterable, 3)
perm = ["".join(ele) for ele in a]

def q_gram_string(DNA_sequence):
    global perm
    #print(perm)
    temp = ''
    for i in perm:
        if i in DNA_sequence:
            temp += '1'
        else:
            temp += '0'
    return temp

def q_gram_sequence(DNA_sequence):
    result = ''
    i = 0
    while i < len(DNA_sequence): #DJORDJE: you should never modify the control variable (i) in the loop like this. You can use while. I'm not sure this is correct. Also, why would you do this in 22-character chunks?
        result += q_gram_string(DNA_sequence[i: i + 20])    #evaluates q=3_gram for consecutive 22 character blocks and appends to results
        i += 20
    #print(result)
    return result

def random_string(alphabet = iterable, length = 4):
    return ''.join(random.choices(iterable, k=length))

def initial_clustering(pool, num_strands):
    clusters = [[i] for i in pool]
    return clusters

def in_clusters(clusters, ele):
    for cluster in clusters:
        if ele in cluster: 
            return True
    return False

def clustering(clusters, global_iter, local_iter, w, l):
    new_clusters = clusters
    #print(clusters)
    for i in range(0, global_iter): #in each global iteration
        if (i==0 and (int(flags['first_iteration_using_index']==1))):    #for the first iteration create partitions using prefix as hash value
            global_partition = prefix_hashing(new_clusters, w, l)
        else:
            global_partition = global_hashing(new_clusters, w, l)
        if(flags.getboolean('print_global_iterations')):
            print("Global iter :", i+1, "#clusters :", len(new_clusters))
        #print(global_partition)
        new_clusters = []
        global_partition = cluster_partitions(global_partition, w, l)
        #print(global_partition)
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
        hash_value = representative[:10]  #DJORDJE: length might need tuning
        if hash_value not in global_partition:
            global_partition[hash_value] = []
        global_partition[hash_value].append(cluster)
    if(flags.getboolean('print_global_partition')):
        count = 1
        for key, cluster in global_partition.items():
            print("partition:", count, "key:", key, "partition size:", len(cluster))
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
        pos = -1
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
            print("partition:", count, "key:", key, "partition size:", len(cluster))
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
        #print("Parition no. :", count, "cluster size :", len(clusters))
        global_partition[key] = local_clustering(clusters, w, l)
        #count += 1
    #print(global_partition)
    return global_partition

def underlying_clustering(pool, actual_num_strands, multiplier):
    clusters = []
    for i in range(0, actual_num_strands):
        temp = []
        for j in range(0, multiplier):
            temp.append(pool[i + actual_num_strands*j])
        clusters.append(temp)
    return clusters

def create_q_gram_dict(pool):
    for ele in pool:
        q_gram_dict[ele] = q_gram_sequence(ele)

def kickWrongEle(clusters, index_length):
    newcluster = []
    wrongEle = []
    
    for ele in clusters:
        compClus = ele
        tempClus = []
        for subele in ele:
            wrongNum = 0
            for compele in compClus:
                dist = compute_distance(subele, compele, "edit")
                if dist > wrong_bases_limit:
                    wrongNum += 1
            if wrongNum > deviation_fraction * len(ele):
                wrongEle.append(subele)
            else:
                tempClus.append(subele)
        if len(tempClus) >= min_cluster_size:
            newcluster.append(tempClus)
        else:
            wrongEle.extend(tempClus)

    return newcluster, wrongEle

def kickSingle(clusters, index_length):
    newcluster = []
    wrongEle = []

    for ele in clusters:
        if len(ele) < min_cluster_size:
            wrongEle.extend(ele)
        else:
            newcluster.append(ele)
    
    return newcluster, wrongEle


def Recluster(cluster, wrongEle, index_length, real_num_strand):
    wrongAgain = []
    for ele in wrongEle:
        correctId = -1
        MINDistance = 100
        for i in range(0, len(cluster)):
            if(len(cluster[i]) < 10):
                representative = random.choice(cluster[i])
                ele_index = ele[0:index_length]
                DNA_index = representative[0:index_length]
                if compute_distance(ele_index, DNA_index, "edit") <= max_edit_limit:
                    dis = compute_distance(ele, representative, "edit")
                    if dis < MINDistance:
                        MINDistance = dis
                        correctId = i
            else:
                continue
        if MINDistance < strand_length*(1-threshold):
            cluster[correctId].append(ele)
        else:
            wrongAgain.append(ele)

    result = [None] * real_num_strand
    if len(cluster) == real_num_strand:
        return cluster, result
    print("wrongAgian length:", len(wrongAgain))

    wrongclusters = initial_clustering(wrongAgain, actual_num_strands)
    sys.stdout = open(os.devnull, 'w')
    wrongclusters = clustering(wrongclusters, global_iter, local_iter, w ,l)
    sys.stdout = sys.__stdout__
    print("Clusters within recycle bin:", len(wrongclusters))
    return cluster + wrongclusters

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
        count += i_count
    accuracy = count/len(underlying_clusters)
    
    return accuracy

def main():
    global gamma
    global actual_num_strands
    
    file_content = read_binary(file_path)

    num_strands, data_length, index_length, data_length_binary = strands_compute(strand_length, len(file_content))

    actual_num_strands = int(len(file_content) / data_length_binary) + 1

    DNA_sequence = encode(file_content, num_strands, data_length, index_length, data_length_binary)

    with open(file_locations['DNA_sequences'], "w") as f:
        for item in DNA_sequence:
            f.write("%s\n" %item)
    
    pool = generate_strand_error(DNA_sequence)
    with open(file_locations['DNA_pool'], "w") as g:
        wr = csv.writer(g)
        wr.writerows(pool)

    create_q_gram_dict(pool)
    pool_time = time.time() - start_time

    underlying_clusters = underlying_clustering(pool, actual_num_strands, 10)

    clusters_paper = initial_clustering(pool, actual_num_strands)

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

    clustering_time = time.time() - pool_time -start_time

    if(flags.getboolean('cleanup') == 1):
        if(flags.getboolean('cleanup_false_positives_or_singles') == 1):
            clusters, wrongEle = kickWrongEle(clusters, index_length)
        else:
            clusters, wrongEle = kickSingle(clusters, index_length)

        print("Clusters after deletion:", len(clusters)) 
        print("Number of strands removed :", len(wrongEle))
        #print(clusters)

        with open(file_locations['cleaned_up_clusters'],"w") as f:
            wr = csv.writer(f)
            wr.writerows(clusters)

        cleanup_time = time.time() - clustering_time - pool_time - start_time
        
        if(flags.getboolean('recycle') == 1):
            clusters = Recluster(clusters, wrongEle, index_length, actual_num_strands)

            print("Number after reclustering :", len(clusters))
        
            with open(file_locations['recycled_clusters'],"w") as f:
                wr = csv.writer(f)
                wr.writerows(clusters)

            recycling_time = time.time() - cleanup_time - clustering_time - pool_time - start_time
    
    print("Accuracy (gamma=", gamma,") :", accuracy(clusters, underlying_clusters))

    accuracy_time = time.time() - clustering_time - pool_time - start_time

    program_execution_time = time.time() - start_time

    if(flags.getboolean('print_time_stamps') == 1):
        print("--- %s seconds for creating pool and q_grams---" % (pool_time))
        print("--- %s seconds for creating underlying and output clusters---" % (clustering_time))
        if(flags.getboolean('cleanup') == 1):
            print("--- %s seconds for cleanup of output clusters---" % (cleanup_time))
            accuracy_time = accuracy_time - cleanup_time
        if(flags.getboolean('recycle') == 1):
            print("--- %s seconds for recycling strands---" % (recycling_time))
            accuracy_time = accuracy_time - recycling_time
        print("--- %s seconds for computing accuracy---" % (accuracy_time))    
        print("--- %s seconds for entire program execution---" % (program_execution_time))    

if __name__=="__main__":
    main()