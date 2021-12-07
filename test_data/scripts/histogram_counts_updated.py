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
f.set_figheight(20)


plt.bar(x, x_counts, log=True, align='center', edgecolor="black")
plt.gca().set_xticks(x)
plt.xlabel("# Substitutions", fontsize=30)
plt.ylabel("# Reads (log scale)",fontsize=30)
plt.title("Reads with number of substitutions", fontsize=50)


for i,j in zip(x,x_counts):
    if j > 99:
        plt.text(i-0.5, j+j*0.3, j, fontsize = 18, rotation=90)
    else:
        plt.text(i-0.25, j+j*0.3, j, fontsize = 16)


plt.show()
