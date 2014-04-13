import numpy as np

l = range(20)
n = 3

chnks = [l[i:i+n] for i in range(0, len(l)-n)]
print chnks