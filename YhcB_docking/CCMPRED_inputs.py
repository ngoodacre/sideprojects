

def make_paired_eggnog_file(wdir,eggnogid1,eggnogid2,prot1,prot2):
    inf1=open(wdir+'\\'+eggnogid1+'_extended_members_'+prot1+'.txt')
    entries1=inf1.read().strip().split('>')[1:]
    inf1.close()
    inf2=open(wdir+'\\'+eggnogid2+'_extended_members_'+prot2+'.txt')
    entries2=inf2.read().strip().split('>')[1:]
    inf2.close()
    outf=open(wdir+'\\'+eggnogid1+'_'+prot1+'_'+eggnogid2+'_'+prot2+'.pairedseqs.fasta','w')
    c=0
    for entry1 in entries1:
        org1=entry1.split('.')[0]
        for entry2 in entries2:
            org2=entry2.split('.')[0]
            if org1==org2:
                c+=1
                header1,seq1=entry1.strip().split('\n')
                header2,seq2=entry2.strip().split('\n')
                outf.write('>'+header1+'__'+header2+'\n'+seq1.strip()+seq2.strip()+'\n')
    outf.close()
    print 'Paired fasta files generated for '+str(c)+' pairs of sequences'

wdir='F:\\yeDUFs\\Jitender'
eggnogid1='COG3105'
prot1='YhcB'
eggnogid2='COG1426'
prot2='RodZ'
eggnogid3='COG1077'
prot3='MreB'

make_paired_eggnog_file(wdir,eggnogid1,eggnogid3,prot1,prot3)
