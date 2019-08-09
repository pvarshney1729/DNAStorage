import random

bases = ['A','C','G','T']
choices = [('sub',4),('in',4),('del',4),('skip',288)]
total = sum(w for c, w in choices)

# http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
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

for iteration in range(2**10):
    strand = []
    for i in range(150):
        strand.append(random.choice(bases))
    # print ''.join(strand)
    for replicate in range(10):
        print add_noise(strand)
