import csv

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

def collect_filenames_main(neg_tags):


def main_loop(align_fasta_filenames,ccmpred_filenames):
    contacts_dict=dict()
    for f,align_fasta_fn in enumerate(align_fasta_filenames):
        fnshort=align_fasta_fn.strip().split('\\')[-1]
        print fnshort
        mat=read_in_raw_ccmpred(ccmpred_filenames[f])
        thresholds_filename='F:\\KINEV\\'+'ccmpreds_cutoffs_skempipair.txt'
        thrd=read_in_thresholds(thresholds_filename)
        thrp=calc_prob_threshold(mat,skempiid,thrd)
        contacts_raw=output_cmmpred_contacts(mat,skempiid,thrp)
        aln_pos_map=get_aln_pos_mapping(align_fasta_fn)
        contacts_ref=map_contacts_uniprot(aln_pos_map,contacts_raw,unipairid)
        contacts_dict[skempiid]=contacts_ref
    return contacts_dict

align_filenames,ccmpred_filenames=collect_filenames_main()
contacts_dict=main_loop(align_filenames,ccmpred_filenames)

