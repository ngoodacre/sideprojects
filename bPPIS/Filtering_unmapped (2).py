all_ppis = open("All_PPIs.txt", "r")
Bact_data = open("Reference_Database.txt","r")
Bact_output = open("Bacterial_PPIs.txt","w")

Bact = {}
Bact_PPI = []

for line in Bact_data:
    line = line.strip()
    Bact[line] = 1

for line in all_ppis:
    line = line.strip()
    line = line.split()

    if line[0] in Bact and line[1] in Bact:
        Bact_PPI.append((line[0], line[1]))

Bact_PPI = set(Bact_PPI)

Bact_output.write('\n'.join('%s\t%s' % x for x in Bact_PPI))

Bact_output.close()