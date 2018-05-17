from math import log

output_contingency = open("contingency_table.txt","w")
Bact_data = open("Bacterial_PPIs.txt","r")
Domain_data = open("C:\Users\Shwetha\Desktop\Protein_Interaction\Step2_Mapping_to_Uniprot\domains.txt","r")
domain = {}
domain_ppi = []
contingency = {}

# making a domain conversion dictionary
for line in Domain_data:
    line = line.strip()
    line = line.split()
    domain[line[0]] = line[1]

# converting ppis into domain pfam ids, filters the ones that are not present in the pfam dictionary
for ppi in Bact_data:
    ppi = ppi.strip()
    ppi = ppi.split()
    if ppi[0] in domain and ppi[1] in domain:
        pairs = (domain[ppi[0]], domain[ppi[1]])
        domain_ppi.append(pairs)
        contingency[pairs] = 0
    #else:
     #   pass
quit()
# making a contingency table(dictionary) with the counts of each domain pair
for i in domain_ppi:
    contingency[i] +=1

# writing the contingency table into output file
for keys in contingency.keys():
    #print '\t', keys[0],
    print >> output_contingency,'\t', keys[0],

#print '\n'
print >> output_contingency, "\n"
i = 0
for keys in contingency.keys():
    #print keys[1], '\t',
    print >> output_contingency, keys[1], '\t',
    for keys, values in contingency.items():
        try:
            #print contingency[keys[0], keys[1]], '\t', # taking the log odds value of the frequencies
            print >> output_contingency, contingency[keys[0], keys[1]], '\t',
        except:
            #print log(0), "\t",
            print >> output_contingency, "0.0\t",
    #print '\n'
    i+=1
    print >> output_contingency, '\n'

output_contingency.close()

print "DONE, check output file. Written %s lines", i