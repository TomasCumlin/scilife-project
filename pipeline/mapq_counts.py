# plots from 0 to 60

import sys
import pysam

batch=sys.argv[2]

samfile = pysam.AlignmentFile(sys.argv[1], 'rb')

x=[]
data=[]

for i in samfile:
    data.append(int(i.mapping_quality))
    if int(i.mapping_quality) not in x:
                x.append(int(i.mapping_quality))

x=sorted(x)
x_counts = []

for i in x:
    x_counts.append(data.count(i))


import matplotlib.pyplot as plt
import numpy as np
cumulative_counts = np.cumsum(x_counts)


f = plt.figure()
f.set_figwidth(30)
f.set_figheight(10)

plt.bar(x, cumulative_counts, log=True, align='center', edgecolor="black")
plt.gca().set_xticks(x)
plt.xlabel("Mapq-value", fontsize=30)
plt.ylabel("# Reads (log scale), Cumulative",fontsize=30)
plt.title(f'# Reads Based on Mapq-value, Batch "{batch}"', fontsize=50, pad=90)
#plt.gca().invert_xaxis()

for i,j in zip(x,cumulative_counts):
    plt.text(i-0.5,j, j, fontsize = 18, rotation=90)
plt.savefig(f'mapq_counts_{batch}.png',bbox_inches='tight')


f = plt.figure()
f.set_figwidth(30)
f.set_figheight(10)


plt.bar(x, cumulative_counts, log=True, align='center', edgecolor="black")
plt.gca().set_xticks(x)
plt.xlabel("Mapq-value", fontsize=30)
plt.ylabel("Cumulative Read Frequency (%) - log scale",fontsize=30)
plt.title(f'Read Frequency Based on Mapq-value, Batch "{batch}"', fontsize=50, pad=90)
#plt.gca().invert_xaxis()

for i,j in zip(x,cumulative_counts):
    plt.text(i-0.5,j,round(j/cumulative_counts[-1]*100,6), fontsize = 18, rotation=90)
plt.savefig(f'mapq_freq_{batch}.png',bbox_inches='tight')
