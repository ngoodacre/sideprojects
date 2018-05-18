import csv
maindir='C:\\Users\\Norman.Goodacre\\Documents\\KINEV'
pdb_chain_uniprot_filename=maindir+'\\'+'pdb_chain_uniprot.csv'
skempi_filename=maindir+'\\'+'SKEMPI_wholedatabase.csv'

def pdb_chain_to_uniprot(pdb_chain_uniprot_filename):
    print "Retrieving PDB id, chain, uniprot id mapping information from: "+pdb_chain_uniprot_filename
    inf=open(pdb_chain_uniprot_filename)
    reader=csv.reader(inf)
    pdb_dict=dict()
    for i,row in enumerate(reader):
        if i<2:
            continue
        pdbid=row[0].strip().upper()
        chain=row[1].strip().upper()
        uniprotid=row[2].strip().upper()
        pdb_dict[pdbid+'_'+chain]=uniprotid
    del reader
    inf.close()
    return pdb_dict

def parse_skempi_pairs(skempi_filename,pdb_dict):
    print "Parsing SKEMPI pairs into Uniprot id pairs"
    inf=open(skempi_filename)
    reader=csv.reader(inf)
    skempi_pairs_pdb=[]
    skempi_pairs=[]
    failed_pdbc=[]
    for i,row in enumerate(reader):
        if i==0:
            continue
        pdbid,chain1,chain2=row[0].strip().split('_')
        pdbid=pdbid.upper()
        chain1=chain1.upper()
        chain2=chain2.upper()
        try:
            prot1=pdb_dict[pdbid+'_'+chain1]
        except KeyError:
            failed_pdbc.append(pdbid+'_'+chain1)
            continue
        try:
            prot2=pdb_dict[pdbid+'_'+chain2]
        except KeyError:
            failed_pdbc.append(pdbid+'_'+chain2)
            continue
        skempi_pairs.append([prot1,prot2])
        skempi_pairs_pdb.append([pdbid+'_'+chain1,pdbid+'_'+chain2])
    del reader
    inf.close()
    return [skempi_pairs,skempi_pairs_pdb,failed_pdbc]

def go_uniprotid_superfamily(query_uniprotid,superfamilytype):
    inf=open(maindir+'\\'+'uniprot_pirsf_iproclass.txt')
    print "Finding "+superfamilytype+" superfamily for protein "+query_uniprotid
    for i,line in enumerate(inf):
        uniprotid,pirid,iproid=line.split('\t')
        if query_uniprotid==uniprotid:
            if superfamilytype=='pirsf':
                if pirid=='':
                    return 'unknown'
                else:
                    return pirid
            if superfamilytype=='iproclass':
                if iproid=='':
                    return 'unknown'
                else:
                    return iproid
    return 'unknown'

def go_uniprotid_superfamily_quick():
    print "Retrieving Uniprot to PIRSF mapping for SKEMPI pairs"
    inf=open(maindir+'\\'+'skempi_prot_pirsf.txt')
    sfdict=dict()
    for i,line in enumerate(inf):
        uniprotid,sfids=line.split('\t')
        uniprotid=uniprotid.strip()
        sfids=sfids.strip().split('; ')
        for s,sfid in enumerate(sfids):
            if sfid=='':
                sfid='unknown'
                sfids[s]=sfid
        sfdict[uniprotid]=sfids
    inf.close()
    return sfdict

def go_superfamily_uniprotids(superfamily,superfamilytype):
    inf=open(maindir+'\\'+'uniprot_pirsf_iproclass.txt')
    prots=[]
    print "Finding uniprot ids for "+superfamily+" "+superfamilytype+" superfamily"
    for i,line in enumerate(inf):
        uniprotid,pirid,iproid=line.split('\t')
        if superfamilytype=='pirsf':
            if pirid==superfamily:
                prots.append(uniprotid)
        if superfamilytype=='iproclass':
            if iproid==superfamily:
                prots.append(uniprotid)
    return prots

def go_skempipairs_superfamily(skempipairs,superfamilytype):
    print "Mapping SKEMPI pairs to "+superfamilytype+" protein superfamilies"
    sfdict=dict()
    for skempipair in skempipairs:
        prot1,prot2=skempipair
        sf1=go_uniprotid_superfamily(prot1,superfamilytype)
        sf2=go_uniprotid_superfamily(prot2,superfamilytype)
        sfdict[prot1]=sf1
        sfdict[prot2]=sf2
    return sfdict

def go_superfamily_superfamilymembers(sfdict,superfamilytype,sftag):
    print "Retrieving superfamily members"
    sfams=set([])
    for prot in sfdict.keys():
        thesepirs=sfdict[prot]
        for thispir in thesepirs:
            if sftag in thispir:
                sfams.add(thispir)
    sfdict2=dict()
    inf=open('F:\\KINEV\\uniprot_pirsf_iproclass.txt')
    for i,line in enumerate(inf):
        uniprotid,pirids,iproids=line.split('\t')
        uniprotid=uniprotid.strip()
        pirids=pirids.strip().split(';')
        pirids=[pirid.strip() for pirid in pirids]
        iproids=iproids.strip().split(';')
        iproids=[iproid.strip() for iproid in iproids]
        if superfamilytype=='pirsf':
            for pirid in pirids:
                if pirid in sfams:
                    try:
                        existing=sfdict2[pirid]
                    except KeyError:
                        existing=[]
                    existing.append(uniprotid)
                    sfdict2[pirid]=existing
        if superfamilytype=='iproclass':
            for iproid in iproids:
                if iproid in sfams:
                    try:
                        existing=sfdict2[iproid]
                    except KeyError:
                        existing=[]
                    existing.append(uniprotid)
                    sfdict2[iproid]=existing
    inf.close()
    return sfdict2

def go_uniprotids_taxids(uniprotids):
    print "Retrieving taxonomic identifiers for "+str(len(uniprotids))+" uniprot ids"
    taxdict=dict()
    inf=open(maindir+'\\'+'uniprot_organisms.txt')
    c=0
    notaxa=[]
    for i,line in enumerate(inf):
        if i-c==1000000:
            print i
            c=i
        try:
            uniprotid,taxid=line.split('\t')
            uniprotid=uniprotid.strip()
            taxid=taxid.strip()
        except ValueError:
            uniprotid=line.strip().split()[0]
            if uniprotid in uniprotids:
                notaxa.append(uniprotid)
            continue
        if uniprotid in uniprotids:
            taxdict[uniprotid]=taxid
    inf.close()
    return [taxdict,notaxa]

def go_uniprotids_taxids_quick():
    print "Collecting uniprot-taxid mapping from file: F:\\KINEV\\skempi_uniprot_sffamily_interactingpairs_allprots_taxidmap.txt"
    inf=open(maindir+'\\'+'skempi_uniprot_sffamily_interactingpairs_allprots_taxidmap.txt')
    taxdict=dict()
    for i,line in enumerate(inf):
        prot,taxid=line.strip().split('\t')
        taxdict[prot]=taxid
    inf.close()
    return taxdict

def go_skempipairs_interacting_superfamily_members(skempi_pairs,sfdict,sfdict2,taxdict):
    print "Retrieving same-taxon potential interacting proteins from superfamilies of SKEMPI pairs"
    idict=dict()
    for skempi_pair in skempi_pairs:
        interactions=[]
        prot1,prot2=skempi_pair
        sf1=sfdict[prot1]
        sf2=sfdict[prot2]
        for sfsub1 in sf1:
            if sfsub1=='unknown':
                continue
            sfprots1=sfdict2[sfsub1]
            for sfprot1 in sfprots1:
                taxid1=taxdict[sfprot1]
                for sfsub2 in sf2:
                    if sfsub2=='unknown':
                        continue
                    sfprots2=sfdict2[sfsub2]
                    for sfprot2 in sfprots2:
                        taxid2=taxdict[sfprot2]
                        if taxid1==taxid2:
                            interactions.append([sfprot1,sfprot2])
        idict['_'.join(skempi_pair)]=interactions
    return idict

def print_skempipdb_skempiuniprot_sffamily_interactingpairs(skempi_pairs_pdb,skempi_pairs,idict):
    print "Generating file with all SKEMPI-Uniprot-PIRSF-Same-species pairs from interacting families"
    outrows=[]
    header=['PDBid1','PDBid2','Uniprotid1','Uniprotid2','ikey','ipairs']
    for p,skempi_pair_pdb in enumerate(skempi_pairs_pdb):
        skempi_pair=skempi_pairs[p]
        skkey='_'.join(skempi_pair)
        fam_pairs=idict[skkey]
        if len(fam_pairs)>=10:
            outrow=[]
            outrow.extend(skempi_pair_pdb)
            outrow.extend(skempi_pair)
            outrow.append(skkey)
            outrow.append('; '.join(['_'.join(ipair) for ipair in fam_pairs]))
            outrows.append(outrow)
    outf=open(maindir+'\\'+'skempi_pairs_uniprot_pairs_sfinteracting_pairs.csv','wb')
    writer=csv.writer(outf)
    writer.writerows(outrows)
    outf.close()

def go_ipairs_fasta():
    print "Generating paired fasta files for all interacting PIRSF family member pairs"
    inf=open(maindir+'\\'+'skempi_pairs_uniprot_pairs_sfinteracting_pairs_allprots.fasta')
    fdict=dict()
    entries=inf.read().strip().split('>')[1:]
    inf.close()
    for entry in entries:
        uniprotid=entry.strip().split('|')[1]
        seq=''.join(entry.strip().split('\n')[1:])
        fdict[uniprotid]=seq
    return fdict

def write_paired_fastafiles_ipairs(skempi_pairs_pdb,skempi_pairs,idict,fdict):
    new_idict=dict()
    for skkey in idict.keys():
        ipairs=idict[skkey]
        if len(ipairs)>=10:
            new_idict[skkey]=ipairs
    print 'writing paired fasta files for '+str(len(new_idict))+' protein pairs from SKEMPI'
    for s,skempi_pair_pdb in enumerate(skempi_pairs_pdb):
        skempi_pair=skempi_pairs[s]
        outfilename=maindir+'\\'+'fasta_files\\'+skempi_pair_pdb[0]+'_'+skempi_pair_pdb[1]+'_'+skempi_pair[0]+'_'+skempi_pair[1]+'.sffamily_interactingpairs.fasta'
        outfile=open(outfilename,'w')
        skkey='_'.join(skempi_pair)
        if not skkey in new_idict.keys():
            continue
        write_paired_fastafile_ipair(skkey,idict,fdict,outfile)
        print 'wrote outfile: '+outfilename

def write_paired_fastafile_ipair(skkey,idict,fdict,outfile):
    ipairs=idict[skkey]
    for ipair in ipairs:
        prot1,prot2=ipair
        seq1=fdict[prot1]
        seq2=fdict[prot2]
        outfile.write('>'+'_'.join(ipair)+'\n'+seq1.strip()+seq2.strip()+'\n')
    outfile.close()

def generated_paired_alignment():
    print "Generating paired alignments in ClustalW format using MUSCLE"
    from Bio.Align.Applications import MuscleCommandline
    muscle_cline = MuscleCommandline(input="F:\\KINEV\\fasta_files\\1FFW_A_1FFW_B_P0AE67_P07363.sffamily_interactingpairs.fasta")
    stdout, stderr = muscle_cline()
    from StringIO import StringIO
    from Bio import AlignIO
    align = AlignIO.read(StringIO(stdout), "fasta")
    print(align)

def read_in_raw_ccmpred(ccmpred_filename):
    print "Reading in CCMPred output file: "+ccmpred_filename
    mat=[]
    inf=open(ccmpred_filename)
    for i,line in enumerate(inf):
        row=line.strip().split()
        row=[float(cell) for cell in row]
        mat.append(row)
    inf.close()
    return mat

def read_in_thresholds(thresholds_filename):
    print "Reading in threshold information from file: "+thresholds_filename
    thrd=dict()
    inf=open(thresholds_filename)
    for i,line in enumerate(inf):
        sl=line.strip().split(', ')
        threshold=float(sl[-1].strip().split()[1].split('%')[0])/100
        sl[-1]=sl[-1].split()[0]
        for pdbchainid in sl:
            thrd[pdbchainid]=threshold
    inf.close()
    return thrd

def calc_prob_threshold(mat,skempiid,thrd):
    print "Calculating coevolving residue pair probability threshold for SKEMPI pair: "+skempiid
    allprobs=[]
    threshold=thrd[skempiid]
    for row in mat:
        allprobs.extend(row)
    allprobs=sorted(allprobs)
    l=len(allprobs)
    ind_thresh=int(threshold*l)
    thrp=allprobs[ind_thresh]
    return thrp

def output_cmmpred_contacts(mat,skempiid,thrp):
    "Calculating contacts for SKEMPI pair: "+skempiid
    contacts=[]
    for row in mat:
        contacts_row=[]
        for cell in row:
            if cell>=thrp:
                contacts_row.append(1)
            else:
                contacts_row.append(0)
        contacts.append(contacts_row)
    return contacts

def get_aln_pos_mapping(aln_filename):
    print "Mapping from alignment to Uniprot position, from file: "+aln_filename
    inf=open(aln_filename)
    entries=inf.read().strip().split('>')[1:]
    posdict=dict()
    for entry in entries:
        posdict_sub=dict()
        unipairid=entry.strip().split('\n')[0].strip()
        alnseq=''.join(entry.strip().split('\n')[1:])
        c=1
        for b,base in enumerate(alnseq):
            if not base=='-':
                posdict_sub[b+1]=c
                c+=1
        posdict[unipairid]=posdict_sub
    return posdict

def map_contacts_uniprot(posdict,contacts,unipairid_filter):
    print "Mapping contacts from alignment to Uniprot position, for Uniprot pair: "+unipairid_filter
    mat_map=[]
    c=0
    for r,row in enumerate(contacts):
        if r-c==10:
            print r
            c=r
        for unipairid in posdict.keys():
            if not unipairid==unipairid_filter:
                continue
            posdict_sub=posdict[unipairid]
            try:
                pos1=posdict_sub[r+1]
            except KeyError:
                continue
            existing_row=[]
            for c,cell in enumerate(row):
                try:
                    pos2=posdict_sub[c+1]
                except KeyError:
                    continue
                if cell==1:
                    existing_row.append(1)
                else:
                    existing_row.append(0)
            mat_map.append(existing_row)
    print r
    return mat_map

##def map_skempiunipair_sfunipair(aln_filename):
##    s_aln_filename=aln_filename.strip().split('.')[0].strip().split('_')
##    skempipair=s_aln_filename[0]+'_'+s_aln_filename[1]+'_'+s_aln_filename[3]
##    skempi_unipairid=s_aln_filename[4]+'_'+s_aln_filename[5]
##    inf=open(aln_filename)
##    entries=inf.read().strip().split('>')[1:]
##    posdict=dict()
##    for entry in entries:
##        unipairid=entry.strip().split('\n')[0].strip()
##        if unipairid==skempi_unipairid:
##            skempialnseq=''.join(entry.strip().split('\n')[1:])
##            break
##    for entry in entries:
##        unipairid=entry.strip().split('\n')[0].strip()
##        if unipairid==skempi_unipairid:
##            continue
##        alnseq=''.join(entry.strip().split('\n')[1:])
##        c1=1
##        c2=1
##        for b1,base1 in enumerate(skempialnseq):
##            base2=alnseq[base1]
##            
##                
##        posdict_sub=dict()

def collect_mutants_skempi(posdict_main,pdbdict):
    inf=open(maindir+'\\'+'SKEMPI_wholedatabase.csv')
    print "Collecting mutant information (position, substitution, affinity)"
    reader=csv.reader(inf)
    mdict=dict()
    for i,row in enumerate(reader):
        if i==0:
            continue
        skempiid=row[0].strip()
        pdbid,chain1,chain2=skempiid.strip().split('_')
        skempi_unipairid=pdbdict[pdbid+'_'+chain1]+'_'+pdbdict[pdbid+'_'+chain2]
        try:
            posdict=posdict_main[skempi_unipairid]
        except KeyError:
            print skempiid
        try:
            existing=mdict[skempiid]
        except KeyError:
            continue
        mut=row[2].strip()
        aff_mut=float(row[6].strip())
        aff_wt=float(row[7].strip())
        mutinfo=[mut,aff_mut,aff_wt]
        existing.append(mutinfo)
        mdict[skempiid]=existing
    inf.close()
    del reader
    return mdict

def intra_vs_inter_protein_contacts(protseq1,protseq2,contacts):
    print "Organizing contacts into intra-protein 1, intra-protein 2, and inter-protein"
    contacts_intra1=[]
    contacts_intra2=[]
    contacts_inter=[]
    l1=len(protseq1)
    l2=len(protseq2)
    for pos1,c1 in enumerate(contacts):
        if pos1+1<=l1:
            intra1_row=[]
            inter_row=[]
        else:
            intra2_row=[]
        for pos2,cell in enumerate(c1):
            if pos2<=pos1:
                continue
            if pos1+1<=l1:
                if pos2+1<=l1:
                    intra1_row.append(cell)
                else:
                    inter_row.append(cell)
            else:
                intra2_row.append(cell)
        if pos1+1<=l1:
            contacts_intra1.append(intra1_row)
            contacts_inter.append(inter_row)
        else:
            contacts_intra2.append(intra2_row)
    return [contacts_intra1,contacts_intra2,contacts_inter]

def collect_filenames_main(neg_tags):
    paired_fasta_dir=maindir+'\\'+'fastainput'
    align_fasta_dir=maindir+'\\'+'MUSCLEout'
    align_aln_dir=maindir+'\\'+'MUSCLE_convert_out'
    ccmpred_dir=maindir+'\\'+'CCMPredout_new'
    paired_fasta_filenames=collect_filenames(paired_fasta_dir,'.fasta',neg_tags)
    align_fasta_filenames=collect_filenames(align_fasta_dir,'.fasta.afa',neg_tags)
    align_aln_filenames=collect_filenames(align_aln_dir,'.fasta.afa.ccmp',neg_tags)
    ccmpred_filenames=collect_filenames(ccmpred_dir,'.fasta.mat',neg_tags)
    return [paired_fasta_filenames,align_fasta_filenames,align_aln_filenames,ccmpred_filenames]

def collect_filenames(wdir,tag,neg_tags):
    import os
    fns=[]
    for fn in os.listdir(wdir):
        flag=False
        for neg_tag in neg_tags:
            if neg_tag in fn:
                flag=True
                break
        if fn.endswith(tag) and not flag:
            fns.append(wdir+'\\'+fn)
    return fns

def main_loop(align_fasta_filenames,ccmpred_filenames,fdict):
    contacts_dict=dict()
    for f,align_fasta_fn in enumerate(align_fasta_filenames):
        fnshort=align_fasta_fn.strip().split('.')[1]
        print align_fasta_fn
        print fnshort
        pdbid=fnshort.split('_')[0][-4:]
        chain1=fnshort.split('_')[1]
        chain2=fnshort.split('_')[3]
        skempiid=pdbid+'_'+chain1+'_'+chain2
        uniprotid1,uniprotid2=fnshort.split('.')[0].split('_')[-2:]
        unipairid='_'.join([uniprotid1,uniprotid2])
        mat=read_in_raw_ccmpred(ccmpred_filenames[f])
        thresholds_filename=maindir+'\\'+'ccmpreds_cutoffs_skempipair.txt'
        thrd=read_in_thresholds(thresholds_filename)
        thrp=calc_prob_threshold(mat,skempiid,thrd)
        contacts_raw=output_cmmpred_contacts(mat,skempiid,thrp)
        aln_pos_map=get_aln_pos_mapping(align_fasta_fn)
        contacts_ref=map_contacts_uniprot(aln_pos_map,contacts_raw,unipairid)
        protseq1=fdict[uniprotid1]
        protseq2=fdict[uniprotid2]
        return contacts_ref
        contacts_intra1,contacts_intra2,contacts_inter=intra_vs_inter_protein_contacts(protseq1,protseq2,contacts_ref)
        cdict=dict()
        cdict['intra1']=contacts_intra1
        cdict['intra2']=contacts_intra2
        cdict['inter']=contacts_inter
        contacts_dict[skempiid]=cdict
    return contacts_dict
             
pdb_dict=pdb_chain_to_uniprot(pdb_chain_uniprot_filename)
skempi_pairs,skempi_pairs_pdb,failed_pdbc=parse_skempi_pairs(skempi_filename,pdb_dict)
##sfdict=go_uniprotid_superfamily_quick()
##sfdict2=go_superfamily_superfamilymembers(sfdict,'pirsf','SF')
##uniprotids=[]
##for sfid in sfdict2.keys():
##    uniprotids.extend(sfdict2[sfid])
##uniprotids=set(uniprotids)
##taxdict=go_uniprotids_taxids_quick()
##idict=go_skempipairs_interacting_superfamily_members(skempi_pairs,sfdict,sfdict2,taxdict)
####print_skempipdb_skempiuniprot_sffamily_interactingpairs(skempi_pairs_pdb,skempi_pairs,idict)
fdict=go_ipairs_fasta()
##write_paired_fastafiles_ipairs(skempi_pairs_pdb,skempi_pairs,idict,fdict)
####generated_paired_alignment()

paired_fasta_filenames,align_fasta_filenames,align_aln_filenames,ccmpred_filenames=collect_filenames_main(['1FFW'])
contacts_dict=main_loop(align_fasta_filenames,ccmpred_filenames,fdict)
protseq1=fdict['P01241']
protseq2=fdict['P10912']
contacts_intra1,contacts_intra2,contacts_inter=intra_vs_inter_protein_contacts(protseq1,protseq2,contacts_dict)

##posdict_main=dict()
##posdict_main[skempiid]=aln_pos_map
##mdict=collect_mutants_skempi(posdict_main,pdb_dict)
