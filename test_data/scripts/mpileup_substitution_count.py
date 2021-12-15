import os

cmd = 'samtools mpileup -C50 -f reference/hg19.with.mt.fasta -r chr1 -l Twist_DNA_ST/pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.bed marked_duplicates.bam'

saving_here = os.popen(cmd).readlines()

file_list=[i.split() for i in saving_here]

def count_substitutions(file_name):

    subst_list=[0]*12
    reference_count=[0]*4

    for i in file_name:
        if i[2]=="A":
            reference_count[0]+=int(i[3])
            subst_list[0]+=int(i[4].count("c"))
            subst_list[0]+=int(i[4].count("C"))
            subst_list[1]+=int(i[4].count("g"))
            subst_list[1]+=int(i[4].count("G"))
            subst_list[2]+=int(i[4].count("t"))
            subst_list[2]+=int(i[4].count("T"))

        if i[2]=="C":
            reference_count[1]+=int(i[3])
            subst_list[3]+=int(i[4].count("a"))
            subst_list[3]+=int(i[4].count("A"))
            subst_list[4]+=int(i[4].count("g"))
            subst_list[4]+=int(i[4].count("G"))
            subst_list[5]+=int(i[4].count("t"))
            subst_list[5]+=int(i[4].count("T"))
        if i[2]=="G":
            reference_count[2]+=int(i[3])
            subst_list[6]+=int(i[4].count("a"))
            subst_list[6]+=int(i[4].count("A"))
            subst_list[7]+=int(i[4].count("c"))
            subst_list[7]+=int(i[4].count("C"))
            subst_list[8]+=int(i[4].count("t"))
            subst_list[8]+=int(i[4].count("T"))
        if i[2]=="T":
            reference_count[3]+=int(i[3])
            subst_list[9]+=int(i[4].count("a"))
            subst_list[9]+=int(i[4].count("A"))
            subst_list[10]+=int(i[4].count("c"))
            subst_list[10]+=int(i[4].count("C"))
            subst_list[11]+=int(i[4].count("g"))
            subst_list[11]+=int(i[4].count("G"))

    reference_count_x = [reference_count[0]]*3+[reference_count[1]]*3+[reference_count[2]]*3+[reference_count[3]]*3

    return subst_list,reference_count, reference_count_x



# obtain nb of subst, nb of reads per reference base, and additional list used for plotting.

subst_data = count_substitutions(file_list)

subst_data_list = [i for i in subst_data[0]]

subst_data_freq = [i/j for i,j in zip(subst_data[0],subst_data[2])]

# Used as x-axis in plots.

subst_string=["A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G",]

# lists of colors for plotting.

colors=["violet"]*3+ ["purple"]*3+["blueviolet"]*3+["blue"]*3

colors_double = ["blue","red","green","darkmagenta","orange","grey","lightgrey","khaki", "violet", "lightgreen","pink","lightblue"]

import matplotlib.pyplot as plt
import numpy as np


# figures:

f = plt.figure()
f.set_figwidth(20)
f.set_figheight(10)

for i,j,k in zip(subst_string,subst_data[0],colors):
    plt.bar(i,j,color=k,edgecolor="black")
    plt.text(i,j+j*0.02,j,ha="center")

plt.xlabel("Type of substitution", fontsize=15)
plt.ylabel("# Substitutions",fontsize=15)
plt.title("Total Substitution Counts", fontsize=20)
plt.savefig('total_substitution_counts.png')
#plt.show()

f = plt.figure()
f.set_figwidth(20)
f.set_figheight(10)

for i,j,k,l in zip(subst_string,subst_data[0],colors,subst_data[2]):
    plt.bar(i,j/l,color=k,edgecolor="black")
    plt.text(i,(j/l)**0.999,str(round(j/l*100,5))+" %",ha="center")

plt.xlabel("Type of substitution", fontsize=15)
plt.ylabel("# Substitutions (%)",fontsize=15)
plt.title("Substitution Frequencies", fontsize=20)
plt.savefig('substitution_frequencies.png')
#plt.show()

f = plt.figure()
f.set_figwidth(20)
f.set_figheight(10)

for i in range(int(len(subst_data_list)/2)):
    a = subst_string[i]+ " or " +subst_string[-(i+1)]
    b = subst_data_list[i]+subst_data_list[-(i+1)]
    plt.bar(a,b,edgecolor="black")
    plt.text(a,b**1.002,b,ha="center")

plt.xlabel("Type of substitutions", fontsize=15)
plt.ylabel("# Substitutions",fontsize=15)
plt.title("Substitution Count", fontsize=20)
plt.savefig("substitution_counts_merged.png")


f = plt.figure()
f.set_figwidth(20)
f.set_figheight(10)

for i in range(int(len(subst_data_list)/2)):
    a = subst_string[i]+ " and " +subst_string[-(i+1)]
    b = subst_data_list[i]
    c = subst_data_list[-(i+1)]
    plt.bar(a,b,edgecolor="black",color=colors_double[i])
    plt.bar(a,c,edgecolor="black", bottom=b,color=colors_double[-(i+1)])
    plt.text(a,b**0.95,("(" + subst_string[i]+ ")"),ha="center")
    plt.text(a,b**0.85,b,ha="center")
    plt.text(a,b+c**0.95,("(" + subst_string[-(i+1)] + ")"),ha="center")
    plt.text(a,b+c**0.85,c,ha="center")
    plt.text(a,b+c**1.01,f"total: {b+c}",ha="center",fontsize = 12)

plt.xlabel("Type of substitutions", fontsize=15)
plt.ylabel("# Substitutions",fontsize=15)
plt.title("Substitution Count", fontsize=20)
plt.savefig("substitution_counts_merged_2")


f = plt.figure()
f.set_figwidth(20)
f.set_figheight(10)

for i in range(int(len(subst_data_freq)/2)):
    a = subst_string[i]+ " and " +subst_string[-(i+1)]
    b = subst_data_freq[i]
    c = subst_data_freq[-(i+1)]
    plt.bar(a,b,edgecolor="black",color=colors_double[i])
    plt.bar(a,c,edgecolor="black", bottom=b,color=colors_double[-(i+1)])
    plt.text(a,b**1.05,("(" + subst_string[i]+ ")"),ha="center")
    plt.text(a,b**1.2,f"{round(b*100,5)} %",ha="center")
    plt.text(a,b+c**1.05,("(" + subst_string[-(i+1)] + ")"),ha="center")
    plt.text(a,b+c**1.2,f"{round(c*100,5)} %",ha="center")
    plt.text(a,b+c**0.994,f"total: {round((b+c)*100,5)} %",ha="center",fontsize = 10)

plt.xlabel("Type of substitutions", fontsize=15)
plt.ylabel("# Substitution (%)",fontsize=15)
plt.title("Substitution Frequencies", fontsize=20)
plt.savefig("substitution_frequencies_merged")
