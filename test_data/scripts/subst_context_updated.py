def read_text_file(plt_file):
    file=open(plt_file, 'r')
    file_read=file.readlines()
    file_list=[i.split() for i in file_read]
    return(file_list)

def context_dependency(file_name):

    subst_list=[0]*12
    reference_count=[0]*4
    bases=["A","a","C","c","G","g","T","t"]
    bases_2=[[bases[0]]*4+[bases[2]]*4+[bases[4]]*4+[bases[6]]*4]
    base_iteration =[[a]*2 for a in range(0,12)]
    base_iteration=[0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11]

    dict_titles_a=[]


    file_name.insert(0,(["sequence start: "]*len(file_name[1])))
    file_name.append(["- sequence end"]*len(file_name[1]))

    for j,i,k in zip(file_name[:-2],file_name[1:-1],file_name[2:]):
        if i[2].upper()==bases[0]:

            reference_count[0]+=int(i[3])

            for y,z in zip(base_iteration[:6],bases[2:]):
                subst_list[y]+=int(i[4].count(z))


        if i[2].upper()==bases[2]:

            reference_count[1]+=int(i[3])

            for y,z in zip(base_iteration[6:12],(bases[:2]+bases[4:])):
                subst_list[y]+=int(i[4].count(z))

        if i[2].upper()==bases[4]:

            reference_count[2]+=int(i[3])

            for y,z in zip(base_iteration[12:18],(bases[:4]+bases[6:])):
                subst_list[y]+=int(i[4].count(z))

        if i[2].upper()==bases[6]:

            reference_count[3]+=int(i[3])

            for y,z in zip(base_iteration[18:],(bases[:-2])):
                subst_list[y]+=int(i[4].count(z))


    reference_count_x = [reference_count[0]]*3+[reference_count[1]]*3+[reference_count[2]]*3+[reference_count[3]]*3

    return subst_list, reference_count, reference_count_x

def plotting_subst(subst_data):

    import matplotlib.pyplot as plt
    import numpy as np


    subst_data_list = [i for i in subst_data[0]]
    subst_data_freq = [i/j for i,j in zip(subst_data[0],subst_data[2])]

    subst_string=["A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G",]
    colors=["violet"]*3+ ["purple"]*3+["blueviolet"]*3+["b"]*3
    colors_double = ["b","r","g","m","y","c","lightcyan","khaki", "violet", "lightgreen","pink","lightblue"]


    f = plt.figure()
    f.set_figwidth(10)

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
        plt.text(a,b+c**1.003,f"total: {b+c}",ha="center",fontsize = 12)

    plt.ylim(0,max(subst_data_list)*2.1)    
    plt.xlabel("Type of substitutions", fontsize=10)
    plt.ylabel("# Substitutions",fontsize=10)
    plt.title("Substitution count", fontsize=15)
    plt.savefig("subst_count.png")
    #plt.show()

    f = plt.figure()
    f.set_figwidth(10)

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
        plt.text(a,b+c**0.994,f"total: {round((b+c)*100,5)}%",ha="center",fontsize = 10)

    plt.ylim(0,max(subst_data_freq)*2.1)
    plt.xlabel("Type of substitutions", fontsize=10)
    plt.ylabel("# Substitution (%)",fontsize=10)
    plt.title("Substitution Frequencies", fontsize=15)
    plt.savefig("subst_freq.png")
   #plt.show()

    f = plt.figure()
    f.set_figwidth(10)

    for i in range(int(len(subst_data_list)/2)):
        a = subst_string[i]+ " or " +subst_string[-(i+1)]
        b = subst_data_list[i]+subst_data_list[-(i+1)]    
        plt.bar(a,b,edgecolor="black")
        plt.text(a,b**1.002,b,ha="center")

    plt.ylim(0,max(subst_data_list)*2.1)
    plt.xlabel("Type of substitutions", fontsize=10)
    plt.ylabel("# Substitutions",fontsize=10)
    plt.title("Substitution count", fontsize=15)
    plt.savefig("subst_count_basic.png")
    #plt.show()

    f = plt.figure()
    f.set_figwidth(10)

    for i,j,k in zip(subst_string,subst_data[0],colors):
        plt.bar(i,j,color=k,edgecolor="black")
        plt.text(i,j+j*0.02,j,ha="center")

    plt.ylim(0,max(subst_data_list)*1.1)
    plt.xlabel("Type of substitution", fontsize=10)
    plt.ylabel("# Substitutions",fontsize=10)
    plt.title("Total substitution counts", fontsize=15)
    plt.savefig("subst_count_total.png")
   #plt.show()

    f = plt.figure()
    f.set_figwidth(10)

    for i,j,k,l in zip(subst_string,subst_data[0],colors,subst_data[2]):
        plt.bar(i,j/l,color=k,edgecolor="black")
        plt.text(i,(j/l)**0.999,str(round(j/l*100,5))+" %",ha="center")

    plt.xlabel("Type of substitution", fontsize=10)
    plt.ylabel("# Substitutions (%)",fontsize=10)
    plt.title("Substitution frequencies", fontsize=15)
    plt.savefig("subst_freq_basic.png")
    #plt.show()

import sys

filen = read_text_file(sys.argv[1])

subst = context_dependency(filen)

plotting_subst(subst)
