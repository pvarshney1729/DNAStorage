import random

strand_len = 250
error_rate = 0.1
sample_cases = 1000

def hamming_distance(s1, s2):
    distance = sum(c1!= c2 for c1,c2 in zip(s1,s2))
    return distance

iterable = ['A', 'G', 'C', 'T']
def random_string(alphabet = iterable, length = 4):
    return ''.join(random.choices(iterable, k=length))

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
    cluster = [original]*5
    rev_cluster = []
    for i in range(0, len(cluster)):
        cluster[i] = add_noise(cluster[i])
        rev_cluster.append(cluster[i][::-1])
    return cluster, rev_cluster

def majority(list):
    return max(set(list), key = list.count)
    
def simple_majority(clu):
    zipped = zip(*clu)
    reconst = ''

    for tup in zipped:
        reconst += max(set(tup), key = tup.count)
    return reconst



def recover_strand(cluster, strand_len):
    ans = ''
    k = 0
    prev_k = -1
    recovered = ''

    for i in range(k, strand_len-1):
        ch_list = refine_majority(cluster, i)
        ch = majority(ch_list)

        '''prev_flag = -1
        if(i>0):
            prev_list = refine_majority(cluster, i-1)
            prev_maj = majority(prev_list)
            prev_flag = 1'''

        for j in range(len(cluster)):
            if len(cluster[j]) == i:
                cluster[j] += ch
                
            if cluster[j][i] != ch:

                
                '''if (prev_flag == -1 or cluster[j][i-1] == prev_maj) and cluster[j][i] == ch2:    #erasure error
                    cluster[j] = cluster[j][:i] + ch + cluster[j][i:]
                elif len(cluster[j]) > i+1:
                    if (prev_flag == -1 or cluster[j][i-1] == prev_maj) and cluster[j][i+1] == ch2: #subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i+1:]
                    if (prev_flag == -1 or cluster[j][i-1] == prev_maj) and cluster[j][i+1] == ch: #insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i+1:]'''
                
                ch2_list = refine_majority(cluster, i+1)
                if len(ch2_list) > 0:
                    ch2 = majority(ch2_list)
                else:
                    ch2 = random.choice(iterable)

                ch3_flag = -1
                if i+2 < strand_len:
                    ch3_flag=1
                    ch3_list = refine_majority(cluster, i+2)
                    if len(ch3_list)>0:
                        ch3 = majority(ch3_list)
                    else:
                        ch3 = random.choice(iterable)

                '''print(i, len(cluster[j]), ch3_flag)
                if cluster[j][i] == ch2 and (ch3_flag==-1 or cluster[j][i+1] == ch3):    #erasure error
                            cluster[j] = cluster[j][:i] + ch + cluster[j][i:]'''
                if len(cluster[j]) > i+2:
                    if cluster[j][i] == ch2 and (ch3_flag==-1 or cluster[j][i+1] == ch3):    #erasure error
                            cluster[j] = cluster[j][:i] + ch + cluster[j][i:]
                    #print(i, len(cluster[j]))
                    elif cluster[j][i+1] == ch2 and (ch3_flag == -1 or cluster[j][1+2] == ch3): #subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i+1:]
                    elif cluster[j][i+1] == ch and cluster[j][i+2] == ch2: #insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i+1:]

                    elif cluster[j][i+1] != ch2:
                        ch4_flag = -1
                        if i+3 < strand_len:
                            ch4_flag=1
                            ch4_list = refine_majority(cluster, i+3)
                            if len(ch4_list)>0:
                                ch4 = majority(ch4_list)
                            else:
                                ch4 = random.choice(iterable)
                        
                        if cluster[j][i] == ch3 and (ch4_flag==-1 or cluster[j][i+1] == ch4):    #erasure error
                                    cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i:]
                        elif len(cluster[j]) > i+3:
                            #print(i, len(cluster[j]))
                            if cluster[j][i+2] == ch3 and (ch4_flag == -1 or cluster[j][1+3] == ch4): #subs
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i+1:]
                            elif cluster[j][i+2] == ch and cluster[j][i+3] == ch2: #insertion error
                                cluster[j] = cluster[j][:i] + cluster[j][i+2:]
                            elif cluster[j][i+1] == ch3 and (ch4_flag == -1 or cluster[j][1+3] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i+1:]
                            elif cluster[j][i+2] == ch2 and cluster[j][i+3] == ch3: 
                                cluster[j] = cluster[j][:i] + ch + cluster[j][i+2:]

                elif len(cluster[j]) == i+2:
                    if cluster[j][i] == ch2:    #erasure error
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i:]
                    #print(i, len(cluster[j]))
                    elif cluster[j][i+1] == ch2: #subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i+1:]
                    elif cluster[j][i+1] == ch: #insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i+1:]
                    else:
                        cluster[j] = cluster[j][:i] + ch               
                        
    
        recovered += ch
        ans = ans[:k] + recovered
        

    last_ch_list = refine_majority(cluster, strand_len - 1)
    
    if len(last_ch_list) > 0:
        last_ch = majority(last_ch_list)
    else:
        last_ch = random.choice(iterable)
    ans += last_ch
    
    return ans

def refine_majority(clu, i):
    ems = []
    for ele in clu:
        if(len(ele) > i):
            ems.append(ele[i])
    if len(ems) == 0:
        ems.append(random.choice(iterable))
    return ems

flag = 0
prob = [0]*strand_len
rev_prob = [0*strand_len]
for i in range(sample_cases):
    original = random_string(iterable, strand_len) 
    cluster, rev_cluster = generate_strand_error(original)
    mj = recover_strand(cluster, strand_len)

    rev_mj = recover_strand(rev_cluster, strand_len)
    rev_rev_mj = rev_mj[::-1]
    mj = mj[0:int(strand_len/2)-1] + rev_rev_mj[int(strand_len/2)-1:strand_len]
    
    for j in range(strand_len):
        if original[j] == mj[j]:
            prob[j] += 1
    #print(i)
    
prob[:] = [x/sample_cases for x in prob]
error = [1-x for x in prob]

print(error)
print("Max error:",max(error))

import matplotlib.pyplot as plt

temp = list(range(1,251))
fig, ax = plt.subplots()
ax.plot(temp, error, label='P = 10%')
legend = ax.legend()

plt.xlabel("Position (1-150)")
plt.ylabel("Probability of incorrect base")
plt.show()
