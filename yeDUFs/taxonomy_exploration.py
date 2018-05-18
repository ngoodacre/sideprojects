taxdir='C:\\Users\\Norman.Goodacre\\taxonomy'
gi_map_filename=taxdir+'\\'+'gi_taxid_nucl.dmp'
names_filename=taxdir+'\\'+'names.dmp'
nodes_filename=taxdir+'\\'+'nodes.dmp'
import re
from collections import defaultdict
import pprint


##################################################################################
########## Commence utility functions                                   ##########
##################################################################################

##### Tree data structure class
def tree():
    return defaultdict(tree)

##### Given a list of taxids, provides you with the Latin names
def map_ncbi_taxids_to_taxnames(names_filename,filterset,short,reverse):
    print 'map_ncbi_taxids_to_taxnames'
    inf=open(names_filename)
    taxdict=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split('\t')
        if not sl[6].strip()=='scientific name':
            continue
        taxid=sl[0].strip()
        if len(filterset)>0:
            if not taxid in filterset:
                continue
        taxname=sl[2].strip()
        if short:
            taxname=' '.join(taxname.split()[:2])
        taxname=taxname.lower()
        if reverse:
            taxdict[taxname]=taxid
        else:
            taxdict[taxid]=taxname
    inf.close()
    return taxdict

########## Initiate taxonomic information ##########
def make_tax_parents(node_filename,filterset):
    print 'make_tax_parents'
    header=['tax_id','parent_tax_id','rank','embl_code','division_id','inherited_div_flag',
            'genetic_code_id','inherited_GC_flag','mitochondiral_genetic_code','inherited_MGC_flag',
            'GenBank_hidden_flag','hidden_subtree_root_flag','comments']
    inf=open(node_filename)
    parents=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split('|')
        sl=[cell.strip() for cell in sl]
        taxid,parent=sl[:2]
        if len(filterset)>0:
            if not taxid in filterset:
                continue
        parents[taxid]=parent
    inf.close()
    return parents

def make_tax_children(node_filename,filterset):
    print 'make_tax_children'
    header=['tax_id','parent_tax_id','rank','embl_code','division_id','inherited_div_flag',
            'genetic_code_id','inherited_GC_flag','mitochondiral_genetic_code','inherited_MGC_flag',
            'GenBank_hidden_flag','hidden_subtree_root_flag','comments']
    inf=open(node_filename)
    children=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split('|')
        sl=[cell.strip() for cell in sl]
        taxid,parent=sl[:2]
        if len(filterset)>0:
            if not taxid in filterset:
                continue
        children[parent]=taxid
    inf.close()
    return children

def make_tax_level(nodes_filename,filterset):
    print 'make_tax_level'
    header=['tax_id','parent_tax_id','rank','embl_code','division_id','inherited_div_flag',
            'genetic_code_id','inherited_GC_flag','mitochondrial_genetic_code','inherited_MGC_flag',
            'GenBank_hidden_flag','hidden_subtree_root_flag','comments']
    inf=open(nodes_filename)
    levels=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split('|')
        sl=[cell.strip() for cell in sl]
        taxid=sl[0]
        rank=sl[2]
        if len(filterset)>0:
            if not taxid in filterset:
                continue
        levels[taxid]=rank
    inf.close()
    return levels









##################################################################################
########## Commence analytical functions                                ##########
##################################################################################

def generate_itol_taxdict():
    inf=open(wdir+'\\ITOL\\'+'iTOL_taxIDs_SpeciesNames.txt')
    itol_taxdict=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split()
        taxid=sl[0]
        itol_name=' '.join(sl[1:])
        itol_taxdict[taxid]=itol_name
    inf.close()
    return itol_taxdict

########## Functions for climbing up taxonomy ##########
##### parents and levels should be unfiltered - i.e. contain mappings for the entire NCBI Taxonomy sheet
def climb_tax_tree_main(taxids,parents,levels,target_level,asdict):
    if asdict:
        target_taxids=dict()
    else:
        target_taxids=[]
    for taxid in taxids:
        try:
            target_taxid=climb_tax_tree(taxid,parents,levels,target_level)
        except KeyError:
            target_taxid='NA'
        if asdict:
            target_taxids[taxid]=target_taxid
        else:
            target_taxids.append(target_taxid)
    return target_taxids

def climb_tax_tree(taxid,parents,levels,target_level):
    level=levels[taxid]
    c=0
    while not level==target_level:
        taxid=parents[taxid]
        level=levels[taxid]
        c+=1
        if c>20:
            taxid='NA'
            break
    return taxid

def get_lineage(taxid,parents,levels,aslevels,reverse):
    level=levels[taxid]
    if aslevels:
        lineage=[level]
    else:
        lineage=[taxid]
    c=0
    while not level=='kingdom':
        taxid=parents[taxid]
        level=levels[taxid]
        if aslevels:
            lineage.append(level)
        else:
            lineage.append(taxid)
        c+=1
        if c>20:
            taxid='NA'
            break
    if not reverse:
        lineage=lineage[::-1]
    return lineage

##### Finds most recent common ancestor - wrapper function #####
def find_MRCA_ext(taxids1,taxids2,parents,levels,mrca_oklevels):
    print "Finding most recent common ancestor between two groups of taxids of lengths "+str(len(taxids1))+" and "+str(len(taxids2))
    all_mrcas=[]
    import operator
    order=operator.itemgetter(-2)
    for taxid1 in taxids1:
        try:
            lvl1=levels[taxid1]
        except KeyError:
            print taxid1
            continue
        mrcas=[]
        for taxid2 in taxids2:
            mrca,depth1,depth2=find_MRCA(taxid1,taxid2,parents,levels)
            mrcas.append([taxid1,taxid2,mrca,depth1,depth2])
        mrcas=sorted(mrcas,key=order)
        keep_mrcas=[]
        for taxid1,taxid2,mrca,depth1,depth2 in mrcas:
            if mrca=='unknown':
                continue
            mrca_lineage_levels=get_lineage(mrca,parents,levels,True,True)
            if len(mrca_oklevels.intersection(set(mrca_lineage_levels)))>0:
                keep_mrcas.append([taxid1,taxid2,mrca,depth1,depth2])
        all_mrcas.append(keep_mrcas)
    return all_mrcas

##### Finds most recent common ancestor - internal function #####
def find_MRCA(taxid1,taxid2,parents,levels):
    reverse=True
    lineage1=get_lineage(taxid1,parents,levels,False,reverse)
    lineage2=get_lineage(taxid2,parents,levels,False,reverse)
    mrca='unknown'
    found=False
    for d1,parent1 in enumerate(lineage1):
        for d2,parent2 in enumerate(lineage2):
            if parent1==parent2:
                mrca=parent1
                found=True
                break
        if found:
            break
    return [mrca,d1,d2]

##### Takes the output of the find_MRCA_ext function together with the iTOL names dict generated by
##### match_iTOL_to_simplified_names and outputs a list of iTOL org names that have been mapped to
def output_iTOL_yedufs(all_mrcas,yeduftaxdict,itoltaxdict):
    print "Selecting final organisms in target taxa, returning yeDUFs for each"
    itol_yedufs=dict()
    for mrcas in all_mrcas:
        if len(mrcas)==0:
            continue
        for mrca in mrcas:
            yeduf_taxid,itol_taxid,mrca_taxid=mrca[:3]
            has_yedufs=yeduftaxdict[0][yeduf_taxid]
            itolname=itoltaxdict[itol_taxid]
            itol_yedufs[itol_taxid]=has_yedufs
    return itol_yedufs

##### Code to write output for all_mrcas
def write_mrca_out(mrcas,outfilename):
    import csv
    outf=open(outfilename,'wb')
    writer=csv.writer(outf)
    outrows=[]
    header=['taxid1','taxid2','taxid mrca','d1','d2']
    outrows.append(header)
    for mrca in mrcas:
        outrows.extend(mrca)
    writer.writerows(outrows)
    outf.close()

filterset=set([])
parents=make_tax_parents(nodes_filename,filterset)
make_tax_children(nodes_filename,filterset)
levels=make_tax_level(nodes_filename,filterset)
infotypes=[]
wdir='C:\\Users\\Norman.Goodacre\\yeDUFs'
from DUF_background import get_taxdict
yeduftaxdict=get_taxdict()
itoltaxdict=generate_itol_taxdict()
##mrca_oklevels=set(['species','species group','species subgroup','subspecies'])
mrca_oklevels=set(['species'])
all_mrcas=find_MRCA_ext(yeduftaxdict[0],itoltaxdict,parents,levels,mrca_oklevels)
itol_yedufs=output_iTOL_yedufs(all_mrcas,yeduftaxdict,itoltaxdict)
