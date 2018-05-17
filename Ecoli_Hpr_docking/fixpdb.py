wdir='C:\\Users\\Norman\\Documents\\PARABON\\example_docking_structures'
infilename=wdir+'\\'+'huLEDGF_vs_hivINTEGRASEmutant9.pdb'
infile=open(infilename)
intext=infile.read().strip().split('ATOM')
infile.close()

newtext=[]
for i,block in enumerate(intext):
    if i==0:
        newtext.append(block)
        continue
    cells=block.strip().split()
    chainid=cells[-1]
    atomtype=cells[1][0]
    if i==len(intext)-1:
        chainid=cells[-2]
        print chainid
    cells.insert(0,'ATOM')
    cells.insert(4,chainid)
    if not i==len(intext)-1:
        cells.pop(-1)
        cells.append(atomtype)
    else:
        cells.pop(-2)
        cells.insert(-1,atomtype)
    try:
        block=cells[0]+' '*(7-len(cells[1]))+cells[1]+'  '+cells[2]+' '*(4-len(cells[2]))+cells[3]+' '+cells[4]+' '*(4-len(cells[5].split('.')[0]))+cells[5]+' '*(8-len(cells[6].split('.')[0]))+cells[6]+' '*(4-len(cells[7].split('.')[0]))+cells[7]+' '*(4-len(cells[8].split('.')[0]))+cells[8]+' '*(3-len(cells[9].split('.')[0]))+cells[9]+' '*(3-len(cells[10].split('.')[0]))+cells[10]+' '*11+cells[11]+'  '
    except IndexError:
        print cells
    if i==len(intext)-1:
        block+='\n'+cells[-1]
    newtext.append(block)
newtext='\n'.join(newtext)

outfilename=infilename.replace('.pdb','_fix.pdb')
outfile=open(outfilename,'w')
outfile.write(newtext)
outfile.close()
