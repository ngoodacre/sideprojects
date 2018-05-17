wdir='F:\\yeDUFs\\Jitender'
aadict=dict(zip(['GLY','ALA','LEU','MET','PHE','TRP','LYS','GLN','GLU','SER','PRO','VAL','ILE','CYS','TYR','HIS','ARG','ASN','ASP','THR'],['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']))
def get_1seq_pdb(pdb_filename,targetchain):
    inf=open(pdb_filename)
    seq=''
    foundnum=set([])
    for line in inf:
        if line.startswith('ATOM'):
            sl=line.strip().split()
            aa_code,chain,resnum=sl[3:6]
            if not chain==targetchain:
                continue
            if resnum in foundnum:
                continue
            else:
                foundnum.add(resnum)
            try:
                aa=aadict[aa_code]
            except KeyError:
                aa='X'
            seq+=aa
    inf.close()
    return seq
            
seq=get_1seq_pdb(wdir+'\\'+'YhcB.pdb','A')
