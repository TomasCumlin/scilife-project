text_file = open("counter_file.txt", "r")

lines = text_file.readlines()


data = [int(i) for i in lines]

x=[]

for i in data:
    if i not in x:
        x.append(i)

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
plt.xlabel("# Cumulative substitutions", fontsize=30)
plt.ylabel("# Reads (log scale)",fontsize=30)
plt.title("Reads with number of substitutions, cumulative", fontsize=50, pad=90)
plt.gca().invert_xaxis()

for i,j in zip(x,cumulative_counts[::-1]):
    plt.text(i+0.5,j+j*0.5, round(j/cumulative_counts[-1]*100,6), fontsize = 18, rotation=90)


plt.show()
