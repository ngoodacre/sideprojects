
def PDB_converter(pdbfilename):
    c=0
    infile=open(pdbfilename)
    outfilename=pdbfilename.split('.pdb')[0]+'_convert.pdb'
    outfile=open(outfilename,'w')
    for line in infile:
        if not line.startswith('ATOM') and not line.startswith('ANISOU'):
            outfile.write(line)
            continue
        sl=line.strip().split()
        chain=sl[4].strip()
        pos=sl[5].strip()
        if chain=='A':
            outfile.write(line)
            continue
        line=line.replace(' '+chain+' ',' A ')
        line=line.replace((4-len(pos))*' '+pos,str(int(pos)+1000))
        outfile.write(line)
    outfile.close()
        
l=PDB_converter('F:\\SNPs@interface\\Docking\\3BLH_PTEFB.pdb')
