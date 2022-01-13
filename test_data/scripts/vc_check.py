def read_text_file(plt_file):
    file=open(plt_file, 'r')
    file_read=file.readlines()
    file_list=[i.split() for i in file_read]
    return(file_list)

import sys

vcf_file = read_text_file(sys.argv[1])


known_vcs = read_text_file(sys.argv[2])


detected_variants = {}
known_vcs_dict = {}

for i in known_vcs[1:-1]:
    if len(i)>1:
        known_vcs_dict[i[1]]=float(i[-1].replace("%",""))
        for j in vcf_file:
            if i[1] in j:
                x = (j[8].split(":"))
                y = (j[9].split(":"))
                if "," in (y[x.index("AF")]):
                    z=y[x.index("AF")].split(",")
                    for index in range(len(z)):
                        detected_variants[j[1]]= round(float(z[index])*100,2)
                else:
                    detected_variants[j[1]]= round((float(y[x.index("AF")])*100),2)

difference = {}
abundance_observed = {}
abundance_expected = {}
un_detected = {}


for i in known_vcs_dict.items():
    if i[0] in detected_variants.keys():
        for j in detected_variants.items():
            if i[0]==j[0]:
                difference[i[0]] = round(abs((i[1]-j[1])/j[1])*100,5)
                abundance_observed[j[0]]=j[1]
                abundance_expected[i[0]]=i[1]
    elif i[0] not in detected_variants.keys():
        un_detected[i[0]]=i[1]

f = open("vc_check_detected.txt", "w")

f.write('Detected variants')
f.write("\n" f'Location Expected abundance Observed Abundance Difference(%)'"\n")

for i,j,k in zip(abundance_expected.items(), abundance_observed.items(), difference.items()):
    f.write(str(i[0])+" "+str(i[1])+" "+str(j[1])+" "+str(k[1])+ "\n")

f.close()


g = open("vc_check_undetected.txt", "w")

g.write('Undetected variants')
g.write("\n" f'Location Expected abundance'"\n")
for i in un_detected.items():
    g.write(str(i[0])+" "+str(i[1])+"\n")

g.close()

h = open("vc_check_summary.txt", "w")

h.write('Summary' "\n")
h.write("# Detected variants: "+str(len(abundance_observed))+"\n"+ "Detected variants (%): " + str(len(abundance_observed)/len(known_vcs[1:-1])*100)+ "\n")
h.write("# Undetected variants: "+ str(len(un_detected))+"\n"+ "Undetected variants (%): " + str(len(un_detected)/len(known_vcs[1:-1])*100)+"\n")
