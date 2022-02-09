# commands: "python <python-file> <vcf-file> <known variants> <"batch name">


# read in files to use
def read_text_file(plt_file):
    file=open(plt_file, 'r')
    file_read=file.readlines()
    file_list=[i.split() for i in file_read]
    return(file_list)


import sys

vcf_file = read_text_file(sys.argv[1])

known_vcs = read_text_file(sys.argv[2])

# saving data in dicts
detected_variants = {}
known_vcs_dict = {}
chromosomes=[]

for i in known_vcs[1:-1]:
    if len(i)>1:
        known_vcs_dict[i[1]]=float(i[-1].replace("%",""))
        chromosomes.append(i[0])
        for j in vcf_file:
            if i[1] in j and i[0] in j:
                x = (j[8].split(":"))
                y = (j[9].split(":"))
                if "," in (y[x.index("AF")]):
                    z=y[x.index("AF")].split(",")
                    for index in range(len(z)):
                        detected_variants[j[1]]= round(float(z[index])*100,2)
                else:
                    detected_variants[j[1]]= round((float(y[x.index("AF")])*100),2)


difference_relative = {}
difference_absolute = {}
abundance_observed = {}
abundance_expected = {}
un_detected = {}
chromosomes_abundant = []
chromosomes_unabundant = []


for i,k in zip(known_vcs_dict.items(),chromosomes):
    if i[0] in detected_variants.keys():
        chromosomes_abundant.append(k)
        for j in detected_variants.items():
            if i[0]==j[0]:
                difference_relative[i[0]] = round(abs((i[1]-j[1])/j[1])*100,5)
                difference_absolute[i[0]]= round((i[1]-j[1]),5)
                abundance_observed[j[0]]=j[1]
                abundance_expected[i[0]]=i[1]
    elif i[0] not in detected_variants.keys():
        un_detected[i[0]]=i[1]
        chromosomes_unabundant.append(k)

batch = sys.argv[3]

f = open(f'vc_check_detected_{batch}.txt', "w")



f.write(f'Chromosome Location Expected(%) Observed(%) Absolute_Difference(%) Relative_Difference(%)'"\n")

for h,i,j,k,l in zip(chromosomes_abundant,abundance_expected.items(), abundance_observed.items(),difference_absolute.items(), difference_relative.items()):
    f.write(str(h)+" "+str(i[0])+" "+str(i[1])+" "+str(j[1])+" "+str(k[1])+" "+str(l[1])+ "\n")

f.close()


g = open(f"vc_check_undetected_{batch}.txt", "w")

g.write(f'Chromosome Location Expected(%)'"\n")
for i,j in zip(chromosomes_unabundant,un_detected.items()):
    g.write(str(i)+" "+str(j[0])+" "+str(j[1])+"\n")

g.close()

h = open(f"vc_check_summary_{batch}.txt", "w")


h.write("# Detected variants: "+str(len(abundance_observed))+"\n"+ "Detected variants (%): " + str(len(abundance_observed)/len(known_vcs[1:-1])*100)+ "\n")
h.write("# Undetected variants: "+ str(len(un_detected))+"\n"+ "Undetected variants (%): " + str(len(un_detected)/len(known_vcs[1:-1])*100)+"\n")

h.close()

#plotting everything
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from sklearn.metrics import r2_score

observed_axis = []
expected_axis = []

for i, j in zip(abundance_observed.items(),abundance_expected.items()):
    if i[0]==j[0]:
        observed_axis.append(i[1])
        expected_axis.append(j[1])

b, m = polyfit(observed_axis,expected_axis, 1)
r2 = r2_score(expected_axis, observed_axis)

reg = []
for i in observed_axis:
    reg.append(b + m * i)

plt.scatter(observed_axis,expected_axis, color="r", alpha=0.5, edgecolors="red")   
plt.plot(observed_axis, reg, '-')
plt.xlabel(f"Observed Abundance (%) of Batch {batch}")
plt.ylabel("Expected Abundance (%)")
plt.title("Expected vs Observed Variant Abundance", fontsize=15)
plt.text(min(observed_axis), max(expected_axis),f"r2={round(r2,3)}",fontsize=12)
plt.savefig(f'vc_check_correlation_{batch}.png',bbox_inches='tight')
