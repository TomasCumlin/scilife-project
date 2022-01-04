import sys
import os

def subst_dicts(pileup_file):

    # creating empty dicts where all mutations will be stored.
    # E.g. "bases_a_c" stores all scenarios where the middle nucleotide is a and mutates into c.

    bases_a_c = {}
    bases_a_g = {}
    bases_a_t = {}
    bases_c_a = {}
    bases_c_g = {}
    bases_c_t = {}
    bases_g_a = {}
    bases_g_c = {}
    bases_g_t = {}
    bases_t_a = {}
    bases_t_c = {}
    bases_t_g = {}

    # freqency dicts
    bases_a_c_freq = {}
    bases_a_g_freq = {}
    bases_a_t_freq = {}
    bases_c_a_freq = {}
    bases_c_g_freq = {}
    bases_c_t_freq = {}
    bases_g_a_freq = {}
    bases_g_c_freq = {}
    bases_g_t_freq = {}
    bases_t_a_freq = {}
    bases_t_c_freq = {}
    bases_t_g_freq = {}


    bases=["A","C","G","T"]
    bases3=["A","C","G","T"]*4
    bases2=["A"]*4+["C"]*4+["G"]*4+["T"]*4

    values=[0]*16

    # Creating all substitutions and possible neghbours, currently only with 0 as the value.

    for i in range(len(bases2)):
        bases_a_c[f'{bases2[i]}A{bases3[i]}->{bases2[i]}C{bases3[i]}'] = values[i]
        bases_a_g[f'{bases2[i]}A{bases3[i]}->{bases2[i]}G{bases3[i]}'] = values[i]
        bases_a_t[f'{bases2[i]}A{bases3[i]}->{bases2[i]}T{bases3[i]}'] = values[i]
        bases_c_a[f'{bases2[i]}C{bases3[i]}->{bases2[i]}A{bases3[i]}'] = values[i]
        bases_c_g[f'{bases2[i]}C{bases3[i]}->{bases2[i]}G{bases3[i]}'] = values[i]
        bases_c_t[f'{bases2[i]}C{bases3[i]}->{bases2[i]}T{bases3[i]}'] = values[i]
        bases_g_a[f'{bases2[i]}G{bases3[i]}->{bases2[i]}A{bases3[i]}'] = values[i]
        bases_g_c[f'{bases2[i]}G{bases3[i]}->{bases2[i]}C{bases3[i]}'] = values[i]
        bases_g_t[f'{bases2[i]}G{bases3[i]}->{bases2[i]}T{bases3[i]}'] = values[i]
        bases_t_a[f'{bases2[i]}T{bases3[i]}->{bases2[i]}A{bases3[i]}'] = values[i]
        bases_t_c[f'{bases2[i]}T{bases3[i]}->{bases2[i]}C{bases3[i]}'] = values[i]
        bases_t_g[f'{bases2[i]}T{bases3[i]}->{bases2[i]}G{bases3[i]}'] = values[i]

    reference_count=[0]*4

    for i,j,k in zip(pileup_file[:-2],pileup_file[1:-1],pileup_file[2:]):
        if j[2].upper()=="A":

            reference_count[0]+=int(j[3])

            for x in j[4]:
                if x.upper()=="C":
                    bases_a_c[f'{i[2].upper()}A{k[2].upper()}->{i[2].upper()}C{k[2].upper()}']+=1
                if x.upper()=="G":
                    bases_a_g[f'{i[2].upper()}A{k[2].upper()}->{i[2].upper()}G{k[2].upper()}']+=1
                if x.upper()=="T":
                    bases_a_t[f'{i[2].upper()}A{k[2].upper()}->{i[2].upper()}T{k[2].upper()}']+=1

        if j[2].upper()=="C":

            reference_count[1]+=int(j[3])

            for x in j[4]:
                if x.upper()=="A":
                    bases_c_a[f'{i[2].upper()}C{k[2].upper()}->{i[2].upper()}A{k[2].upper()}']+=1
                if x.upper()=="G":
                    bases_c_g[f'{i[2].upper()}C{k[2].upper()}->{i[2].upper()}G{k[2].upper()}']+=1
                if x.upper()=="T":
                    bases_c_t[f'{i[2].upper()}C{k[2].upper()}->{i[2].upper()}T{k[2].upper()}']+=1


        if j[2].upper()=="G":

            reference_count[2]+=int(j[3])

            for x in j[4]:
                if x.upper()=="A":
                    bases_g_a[f'{i[2].upper()}G{k[2].upper()}->{i[2].upper()}A{k[2].upper()}']+=1
                if x.upper()=="C":
                    bases_g_c[f'{i[2].upper()}G{k[2].upper()}->{i[2].upper()}C{k[2].upper()}']+=1
                if x.upper()=="T":
                    bases_g_t[f'{i[2].upper()}G{k[2].upper()}->{i[2].upper()}T{k[2].upper()}']+=1


        if j[2].upper()=="T":

            reference_count[3]+=int(j[3])

            for x in j[4]:
                if x.upper()=="A":
                    bases_t_a[f'{i[2].upper()}T{k[2].upper()}->{i[2].upper()}A{k[2].upper()}']+=1
                if x.upper()=="C":
                    bases_t_c[f'{i[2].upper()}T{k[2].upper()}->{i[2].upper()}C{k[2].upper()}']+=1
                if x.upper()=="G":
                    bases_t_g[f'{i[2].upper()}T{k[2].upper()}->{i[2].upper()}G{k[2].upper()}']+=1


    # converting counts to frequences and store them in other dicts.

    for i,j,k in zip(bases_a_c.items(),bases_a_g.items(),bases_a_t.items()):
        bases_a_c_freq[i[0]]=round(i[1]/reference_count[0]*100,14)
        bases_a_g_freq[j[0]]=round(j[1]/reference_count[0]*100,14)
        bases_a_t_freq[k[0]]=round(k[1]/reference_count[0]*100,14)

    for i,j,k in zip(bases_c_a.items(),bases_c_g.items(),bases_c_t.items()):
        bases_c_a_freq[i[0]]=round(i[1]/reference_count[1]*100,14)
        bases_c_g_freq[j[0]]=round(j[1]/reference_count[1]*100,14)
        bases_c_t_freq[k[0]]=round(k[1]/reference_count[1]*100,14)

    for i,j,k in zip(bases_g_a.items(),bases_g_c.items(),bases_g_t.items()):
        bases_g_a_freq[i[0]]=round(i[1]/reference_count[2]*100,14)
        bases_g_c_freq[j[0]]=round(j[1]/reference_count[2]*100,14)
        bases_g_t_freq[k[0]]=round(k[1]/reference_count[2]*100,14)

    for i,j,k in zip(bases_t_a.items(),bases_t_c.items(),bases_t_g.items()):
        bases_t_a_freq[i[0]]=round(i[1]/reference_count[3]*100,14)
        bases_t_c_freq[j[0]]=round(j[1]/reference_count[3]*100,14)
        bases_t_g_freq[k[0]]=round(k[1]/reference_count[3]*100,14)

    subst_counts = [bases_a_c, bases_a_g, bases_a_t, bases_c_a, bases_c_g, bases_c_t, bases_g_a, bases_g_c, bases_g_t, bases_t_a, bases_t_c, bases_t_g]
    subst_freq = [bases_a_c_freq, bases_a_g_freq, bases_a_t_freq, bases_c_a_freq, bases_c_g_freq, bases_c_t_freq, bases_g_a_freq, bases_g_c_freq, bases_g_t_freq, bases_t_a_freq, bases_t_c_freq, bases_t_g_freq]

    transition = [bases_a_g, bases_c_t, bases_g_a, bases_t_c]

    transversion = [bases_a_c, bases_a_t, bases_c_a, bases_c_g, bases_g_c, bases_g_t, bases_t_a, bases_t_g]

    transition_freq = [bases_a_g_freq, bases_c_t_freq, bases_g_a_freq, bases_t_c_freq]

    transversion_freq = [bases_a_c_freq, bases_a_t_freq, bases_c_a_freq, bases_c_g_freq, bases_g_c_freq, bases_g_t_freq, bases_t_a_freq, bases_t_g_freq]    

    return subst_counts, subst_freq, transition, transversion, transition_freq, transversion_freq


def plotting_subst_neighb(bases, freq=False):

    import matplotlib.pyplot as plt

    subst_string=["A->C","A->G","A->T","C->A","C->G","C->T","G->A","G->C","G->T","T->A","T->C","T->G"]

    for x in range(int(len(bases)/2)):

        f = plt.figure()
        f.set_figwidth(30)


        for i,j in zip(bases[x].items(),list(reversed(sorted(bases[-(x+1)].items())))):
            a = j[0]+ "\n and \n" + i[0]
            b = i[1]
            c = j[1]
            plt.bar(a,b,color="lightblue",edgecolor="black")
            plt.bar(a,c, bottom=i[1],color="pink",edgecolor="black")
            if freq is True:
                plt.text(a,b+c**0.99,f"{round((b+c),3)} %",ha="center",fontsize = 12)
            else:
                plt.text(a,b+c**1.01,f"total: {b+c}",ha="center",fontsize = 12)


        plt.ylim(0,(max(bases[x].values())+max(bases[-(x+1)].values()))*1.3)
        plt.xlabel("Type of substitutions", fontsize=20)

        if freq is True:
            plt.ylabel("Substitutions(%)",fontsize=20)
            plt.title(f'Substitution frequencies for {subst_string[x]} and {subst_string[-(x+1)]}', fontsize=25)
            plt.savefig(f'subst_neighbour_freq_{x}.png',bbox_inches='tight')
        else:
            plt.ylabel("# Substitutions",fontsize=20)
            plt.title(f'Substitution count for {subst_string[x]} and {subst_string[-(x+1)]}', fontsize=25)
            plt.savefig(f'subst_neighbour_{x}.png',bbox_inches='tight')


def plotting_subst(subst_data):

    import matplotlib.pyplot as plt
    import numpy as np


    subst_data_list = [sum(i.values()) for i in subst_data[0]]
    subst_data_freq = [sum(i.values()) for i in subst_data[1]]

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
    plt.savefig(f'subst_count',bbox_inches='tight')

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
    plt.savefig(f'subst_freq',bbox_inches='tight')

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
    plt.savefig(f'subst_count_plain',bbox_inches='tight')

def trans_ratio(dicts_data_1, dicts_data_2,dicts_data_3,dicts_data_4):

    transi = 0
    transv = 0
    transi_freq = 0
    transv_freq = 0

    for i in dicts_data_1:
        transi += sum(i.values())

    for i in dicts_data_2:
        transv += sum(i.values())

    for i in dicts_data_3:
        transi_freq += sum(i.values())

    for i in dicts_data_4:
        transv_freq += sum(i.values())

    ratio_counts = transi/transv
    ratio_counts_freq = transi_freq/transv_freq

    return ratio_counts, ratio_counts_freq



cmd = 'samtools mpileup -C50 -f reference/hg19.with.mt.fasta -l Twist_DNA_ST/pool1_pool2_nochr_3c.sort.merged.padded20.hg19.210311.met.annotated.bed marked_duplicates.bam'

pileup_file = os.popen(cmd).readlines()

file_list=[i.split() for i in pileup_file]

file_list.insert(0,(["C"]*len(file_list[1])))
file_list.append(["C"]*len(file_list[1]))

data_dicts = subst_dicts(file_list)

plotting_subst(data_dicts)

plotting_subst_neighb(data_dicts[0])
plotting_subst_neighb(data_dicts[1], freq=True)


transition = trans_ratio(data_dicts[2], data_dicts[3], data_dicts[4], data_dicts[5])
print(f'The transition/transversion ratio: {round(transition[0],4)} (counts) or {round(transition[1],4)} (frequency)')
