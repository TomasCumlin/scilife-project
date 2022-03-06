import sys
import pysam

batch=sys.argv[2]

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

counter = 0
x=[]
data=[]

for read in samfile:
    if "MD" in str(read.tags[1]):
        counter += str(read.tags[1]).count("C")
        counter += str(read.tags[1]).count("c")
        counter += str(read.tags[1]).count("G")
        counter += str(read.tags[1]).count("g")
        counter += str(read.tags[1]).count("A")
        counter += str(read.tags[1]).count("a")
        counter += str(read.tags[1]).count("T")
        counter += str(read.tags[1]).count("t")
        data.append(counter)
        if counter not in x:
            x.append(counter)
        counter = 0

samfile.close()

x=sorted(x)
x_counts = []

for i in x:
    x_counts.append(data.count(i))

import matplotlib.pyplot as plt
import numpy as np

f = plt.figure()
f.set_figwidth(30)
f.set_figheight(10)


cumulative_counts = np.cumsum(x_counts[::-1])

plt.bar(x, cumulative_counts[::-1], log=True, align='center', edgecolor="black")
plt.gca().set_xticks(x)
plt.xlabel("# Substitutions", fontsize=30)
plt.ylabel("Cumulative Reads (%) - Log Scale",fontsize=30)
plt.title(f'Read Frequency Based on # Substitutions, Batch "{batch}"', fontsize=50, pad=90)
plt.gca().invert_xaxis()

for i,j in zip(x,cumulative_counts[::-1]):
    plt.text(i+0.5,j+j*0.5, round(j/cumulative_counts[-1]*100,6), fontsize = 18, rotation=90)


plt.savefig(f"subst_cumulative_{batch}.png",bbox_inches='tight')
