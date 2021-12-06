text_file = open("counter_file.txt", "r")

lines = text_file.readlines()


data = [int(i) for i in lines]

x=[] 

for i in data:
    if i not in x:
        x.append(i)
        
x=sorted(x)        

import matplotlib.pyplot as plt
import numpy as np

        
f = plt.figure()
f.set_figwidth(20)
bin_size = len(x)

plt.hist(data, bins = bin_size, density=True, ec="black", log=True)
#plt.xlim(2,10)
plt.xlabel("# Substitions")
plt.ylabel("# Reads (log scale)")
plt.title("Reads with number of substituions")
plt.xticks(x)