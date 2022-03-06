
# import bed-file

def read_file(file):
    files = open(file,"r")
    files = files.read().split('\n')
    list1 = [i.split() for i in files]
    return list1

# plot coverage of the data

def plot_coverage(data_org,batch):

    colors = {
            "chr1": "r",
            "chr2": "b",
            "chr3": "g",
            "chr4": "y",
            "chr5": "m",
            "chr6": "c",
            "chr7": "lightcoral",
            "chr8": "deepskyblue",
            "chr9": "forestgreen",
            "chr10": "gold",
            "chr11": "hotpink",
            "chr12": "darkcyan",
            "chr13": "tomato",
            "chr14": "cornflowerblue",
            "chr15": "lime",
            "chr16": "chocolate",
            "chr17": "plum",
            "chr18": "gray",
            "chr19": "indianred",
            "chr20": "deepskyblue",
            "chr21": "darkseagreen",
            "chr22": "khaki",
            "chrX": "lightpink",
            "chrY": "lightblue"
        }


    import matplotlib.pyplot as plt

    counter=0
    maximum=0


    f = plt.figure()
    f.set_figwidth(20)


    for i in data_org:
        if len(i)==5:
            if "CNV" not in i[3]:
                plt.bar(counter,float(i[-1]),color=colors[i[0]])
                counter+=1
                if maximum < float(i[-1]):
                    maximum = float(i[-1])

    intervall=round(counter/23,0)
    counter2=0

    for i in colors.items():
        t = plt.text(counter2,maximum*0.9,f'{i[0]}')
        t.set_bbox(dict(facecolor=i[1]))
        counter2+=intervall

    plt.ylim(0,maximum)
    plt.xticks([])
    plt.xlabel("Position", fontsize=11,labelpad=10)
    plt.ylabel("Coverage (Counts)",fontsize=10)
    plt.title(f"Coverage of Batch {batch}", fontsize=15)
    plt.savefig(f'coverage_{batch}.png', bbox_inches='tight')


def difference_relative(data1,data2,batch):

    import numpy as np

    difference = []
    difference_count=[]

    for i,j in zip(data1, data2):
        if len(i) == 5 and len(j) == 5:
            if "CNV" not in i[3]:
                if i[1]==j[1] and i[2]==j[2]:
                    if (float(i[-1])-float(j[-1]))==0:
                        difference.append([i[0],i[1],i[2],i[3],i[-1],j[-1],round((float(i[-1])-float(j[-1])),6),0])
                        difference_count.append([i[0],i[1],i[2],i[3],i[-1],j[-1],round((float(i[-1])-float(j[-1])),6),0])
                    else:
                        difference.append([i[0],i[1],i[2],i[3],i[-1],j[-1],round((float(i[-1])-float(j[-1])),6),round((float(i[-1])-float(j[-1]))/(float(i[-1]))*100,6)])
                        difference_count.append([i[0],i[1],i[2],i[3],i[-1],j[-1],round((float(i[-1])-float(j[-1])),6),round((float(i[-1])-float(j[-1]))/(float(i[-1]))*100,6)])


    length = len(difference)
    for i in range(0, length):
        for j in range(0, length-i-1):
            if (difference[j][-1] < difference[j + 1][-1]):
                temp = difference[j]
                difference[j]= difference[j + 1]
                difference[j + 1]= temp

    length = len(difference_count)
    for i in range(0, length):
        for j in range(0, length-i-1):
            if (difference_count[j][-2] < difference_count[j + 1][-2]):
                temp = difference_count[j]
                difference_count[j]= difference_count[j + 1]
                difference_count[j + 1]= temp

    f = open(f"bed_differences_{batch}.txt", "w")
    f.write(f'Chromosome Start End Name Initial_Abundance New_Abundance Absolute_Difference Relative_Difference(%)\n')

    for i in difference:
        f.write(f'{i[0]} {i[1]} {i[2]} {i[3]} {i[4]} {i[5]} {i[6]} {i[7]}\n')

    f.close()

    return difference, difference_count

# plotting the regions with highest absolute difference from untrimmed data

def cov_diff_top(the_difference,batch):

    import matplotlib.pyplot as plt


    colors = {
            "chr1": "r",
            "chr2": "b",
            "chr3": "g",
            "chr4": "y",
            "chr5": "m",
            "chr6": "c",
            "chr7": "lightcoral",
            "chr8": "deepskyblue",
            "chr9": "forestgreen",
            "chr10": "gold",
            "chr11": "hotpink",
            "chr12": "darkcyan",
            "chr13": "tomato",
            "chr14": "cornflowerblue",
            "chr15": "lime",
            "chr16": "chocolate",
            "chr17": "plum",
            "chr18": "gray",
            "chr19": "indianred",
            "chr20": "deepskyblue",
            "chr21": "darkseagreen",
            "chr22": "khaki",
            "chrX": "lightpink",
            "chrY": "lightblue"
        }


    counter=0
    maximum=0


    f = plt.figure()
    f.set_figwidth(14)
    f.set_figheight(14)

    for i in the_difference[1]:
        plt.barh(counter,float(i[-2]), color=colors[i[0]], edgecolor="0.1")
        plt.text((float(i[-2]))+1,counter,f"{float(i[-2])}",weight='bold')
        plt.text(-70,counter-0.1,i[0],ha="right", fontsize=10)
        plt.text(-70, counter+0.2,f'{i[1]} - {i[2]}',ha="right", fontsize=8)
        counter+=1

        if maximum < float(i[-2]):
            maximum = float(i[-2])
        if counter == 15:
            break

    plt.xlim(0,maximum*1.3)
    plt.gca().invert_yaxis()
    plt.yticks([])
    plt.xlabel("Absolute Difference", fontsize=14,labelpad=10)
    plt.ylabel("Positions",fontsize=14,labelpad=130)
    plt.title(f'Coverage Difference between untrimmed data and Batch "{batch}"', fontsize=15)
    plt.savefig(f'highest_cov_diff_{batch}.png', bbox_inches='tight')

def plot_all_diff(the_difference,batch):

    import matplotlib.pyplot as plt

    counter=0
    maximum=0


    f = plt.figure()
    f.set_figwidth(14)
    f.set_figheight(14)

    for i in the_difference[1]:
        plt.barh(counter,float(i[-2]), color="slateblue")
        counter+=1


    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.yticks([])
    plt.xlabel("Absolute Difference (Logarithmic Scale)", fontsize=14)
    plt.ylabel("Rank of position",fontsize=14)
    plt.title(f'Coverage Difference between untrimmed data and Batch "{batch}"', fontsize=15)
    plt.savefig(f'all_cov_dif_{batch}.png', bbox_inches='tight')


import sys
import os

batch_name=sys.argv[3]

data_untrimmed = read_file(sys.argv[1])
data_comp = read_file(sys.argv[2])

#plot_coverage(data_untrimmed,batch_name)

diff=difference_relative(data_untrimmed,data_comp,batch_name)

cov_diff_top(diff,batch_name)
plot_all_diff(diff,batch_name)
