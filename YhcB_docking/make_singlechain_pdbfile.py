

def monomerize_polymer(pdbfilename):
    import re
    beginfinder=re.compile('ATOM\s+|ANISOU\s+')
    chain_and_resfinder=re.compile('\s+[A-Z]{3}\s+[A-Z]{1}\s+[0-9]+')
    inf=open(pdbfilename)
    outf=open(pdbfilename.replace('.pdb','.c1.pdb'),'w')
    atomline=False
    foundchains=set([])
    atom_bnum=0
    res_bnum=0
    for line in inf:
        if line=='TER\n':
            continue
        try:
            begin=beginfinder.findall(line)[0]
            if line.startswith(begin):
                atomline=True
            else:
                atomline=False
        except IndexError:
            atomline=False
        if atomline:
            sl=line.split()
            chain=sl[4]
            try:
                chain_and_res=chain_and_resfinder.findall(line)[0]
            except IndexError:
                print line
            if not chain=='A':
                if not chain in foundchains:
                    atom_bnum+=atom_num
                    res_bnum+=res_num
                    foundchains.add(chain)
            atom_num=int(sl[1])
            res_num=int(sl[5])
            atom_newnum=atom_bnum+atom_num
            res_newnum=res_bnum+res_num
            line=line.replace(begin+str(atom_num),begin+str(atom_newnum))
            chain_and_res_new=chain_and_res.replace(' '+chain+' ',' A ')
            chain_and_res_new=chain_and_res_new.replace(str(res_num)+' ',str(res_newnum))
            line=line.replace(chain_and_res,chain_and_res_new)
            outf.write(line)
        else:
            if line=='END\n':
                outf.write('TER\nEND\n')
                outf.close()
            else:
                outf.write(line)
    inf.close()

pdbfilename='F:\\yeDUFS\\Jitender\\YhcB.pdb'
monomerize_polymer(pdbfilename)
