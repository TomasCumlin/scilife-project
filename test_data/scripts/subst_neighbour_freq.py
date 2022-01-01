def read_text_file(plt_file):
    file=open(plt_file, 'r')
    file_read=file.readlines()
    file_list=[i.split() for i in file_read]

    # adding an empty row to beginning and the end.
    file_list.insert(0,(["C"]*len(file_list[1])))
    file_list.append(["C"]*len(file_list[1]))

    return(file_list)


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

    return subst_counts, subst_freq


def plotting_subst(bases,freq=False):

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

import sys

file_name = read_text_file(sys.argv[1])
the_data = subst_dicts(file_name)
plotting_subst(the_data[0])
plotting_subst(the_data[1],freq=True)
