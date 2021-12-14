#First code draft on how to plot the pileup-data.

def read_text_file(plt_file, make_list="y"):
    file=open(plt_file, 'r')
    file_read=file.readlines()

    if make_list=="y":
        file_list=[i.split() for i in file_read]
    else:
        file_list=file_read

    return(file_list)

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


filen = read_text_file('mpileup_testing.pileup')

subst_data = count_substitutions(filen)

subst_string=["A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G",]

colors=["violet"]*3+ ["purple"]*3+["blueviolet"]*3+["blue"]*3

import matplotlib.pyplot as plt
import numpy as np

f = plt.figure()
f.set_figwidth(10)

for i,j,k in zip(subst_string,subst_data[0],colors):
    plt.bar(i,j,color=k,edgecolor="black")
    plt.text(i,j+j*0.02,j,ha="center")

plt.xlabel("Type of substitution", fontsize=10)
plt.ylabel("# Substitutions",fontsize=10)
plt.title("Total substitution counts", fontsize=15)
plt.show()

f = plt.figure()
f.set_figwidth(10)

for i,j,k,l in zip(subst_string,subst_data[0],colors,subst_data[2]):
    plt.bar(i,j/l,color=k,edgecolor="black")
    plt.text(i,(j/l)**0.999,str(round(j/l*100,5))+" %",ha="center")

plt.xlabel("Type of substitution", fontsize=10)
plt.ylabel("# Substitutions (%)",fontsize=10)
plt.title("Substitution frequencies", fontsize=15)
plt.show()
