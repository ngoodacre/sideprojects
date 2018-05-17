
def COG_sequences_byspecies_paired(pairs,sdict1,sdict2,paired_fasta_outfn):
    outf=open(paired_fasta_outfn,'w')
    for org in pairs.keys():
        subpairs=pairs[org]
        for subpair in subpairs:
            name1,name2=subpair
            seq1=sdict1[name1]
            seq2=sdict2[name2]
            outf.write('>'+name1+'::'+name2+'('+org+')'+'\n'+seq1+seq2+'\n')
    outf.close()

def COG_sdict(COG_seqs_fn):
    inf=open(COG_seqs_fn)
    entries=inf.read().strip().split('>')[1:]
    sdict=dict(zip([entry.strip().split('\n')[0].strip() for entry in entries],
                   [entry.strip().split('\n')[1].strip() for entry in entries]))
    inf.close()
    return sdict

def COG_members_byspecies_paired(COG1_members_fn,COG2_members_fn):
    inf1=open(COG1_members_fn)
    inf2=open(COG2_members_fn)
    entries1=inf1.read().strip().split('\n')
    entries2=inf2.read().strip().split('\n')
    inf1.close()
    inf2.close()
    names1=[entry1.strip().split('\t')[3].strip()+'.'+entry1.strip().split('\t')[0].strip() for entry1 in entries1]
    names2=[entry2.strip().split('\t')[3].strip()+'.'+entry2.strip().split('\t')[0].strip() for entry2 in entries2]
    orgs1=[entry1.strip().split('\t')[2] for entry1 in entries1]
    orgs2=[entry2.strip().split('\t')[2] for entry2 in entries2]
    mdict1=COG_members_byspecies(orgs1,names1)
    mdict2=COG_members_byspecies(orgs2,names2)
    pairs=dict()
    for org in mdict1.keys():
        names1=mdict1[org]
        try:
            names2=mdict2[org]
        except KeyError:
            continue
        subpairs=[]
        for name1 in names1:
            for name2 in names2:
                subpairs.append([name1,name2])
        pairs[org]=subpairs
    return pairs

##### Auxiliary function of the above
def COG_members_byspecies(orgs,names):
    mdict=dict()
    for o,org in enumerate(orgs):
        name=names[o]
        try:
            existing=mdict[org]
        except KeyError:
            existing=[]
        existing.append(name)
        mdict[org]=existing
    return mdict

wdir='C:\\Users\\Norman.Goodacre\\Documents\\Ecoli_Hpr_docking'       
COG1='COG1925'
COG2='COG0563'
COG1_members_fn=wdir+'\\'+COG1+'_extended_members2.txt'
COG2_members_fn=wdir+'\\'+COG2+'_extended_members2.txt'
COG1_seqs_fn=wdir+'\\'+COG1+'_extended_members.txt'
COG2_seqs_fn=wdir+'\\'+COG2+'_extended_members.txt'
pairs=COG_members_byspecies_paired(COG1_members_fn,COG2_members_fn)
sdict1=COG_sdict(COG1_seqs_fn)
sdict2=COG_sdict(COG2_seqs_fn)
paired_fasta_outfn=wdir+'\\'+COG1+'_'+COG2+'_'+'byspecies_paired.fasta'
COG_sequences_byspecies_paired(pairs,sdict1,sdict2,paired_fasta_outfn)
