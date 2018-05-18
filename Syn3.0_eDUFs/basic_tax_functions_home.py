taxdir='C:\\Users\\Norman\\Documents\\VCU_research\\taxonomy'
gi_map_filename=taxdir+'\\'+'gi_taxid_nucl.dmp'
names_filename=taxdir+'\\'+'taxdump'+'\\'+'names.dmp'
nodes_filename=taxdir+'\\'+'taxdump'+'\\'+'nodes.dmp'
import re
from collections import defaultdict
import pprint

##### Tree data structure class
def tree():
    return defaultdict(tree)

##### Given a raw RefSeq viral genomic file, removes phage and formats so that DNA sequence
##### is on one line only
def trim_refseq_rawfile(refseq_raw_fastafn,refseq_trim_fastafn):
    inf=open(refseq_raw_fastafn)
    entries=inf.read().strip().split('>gi')[1:]
    inf.close()
    outf=open(refseq_trim_fastafn,'w')
    for entry in entries:
        lines=entry.strip().split('\n')
        header=lines[0]
        if ' phage ' in header.lower():
            continue
        seq=''.join([seqline.strip() for seqline in lines[1:]])
        outf.write('>gi'+header+'\n'+seq+'\n')
    outf.close()

##### Given a fasta file, gives you the gis
def get_gis_fastafile(fasta_filename):
    gis=set([])
    inf=open(fasta_filename)
    for i,line in enumerate(inf):
        if line.startswith('>gi'):
            gi=line.split('|')[1]
            gis.add(gi)
    inf.close()
    return gis

##### Given a list of gis, provides you with the NCBI taxids
def map_gis_to_taxids(gis):
    c=0
    mapdict=dict()
    inf=open(gi_map_filename)
    for i,line in enumerate(inf):
        if i-c==1000000:
            print i
            c=i
        gi,taxid=line.strip().split()
        if gi in gis:
            mapdict[gi]=taxid
    inf.close()
    return mapdict

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
            'genetic_code_id','inherited_GC_flag','mitochondiral_genetic_code','inherited_MGC_flag',
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

def get_lineage(taxid,parents,levels,reverse):
    lineage=[taxid]
    level=levels[taxid]
    c=0
    while not level=='kingdom':
        taxid=parents[taxid]
        level=levels[taxid]
        lineage.append(taxid)
        c+=1
        if c>20:
            taxid='NA'
            break
    if not reverse:
        lineage=lineage[::-1]
    return lineage

##### Finds most recent common ancestor - wrapper function #####
def find_MRCA_ext(taxids1,taxids2,parents,levels):
    print "Finding most recent common ancestor between two groups of taxids of lengths "+str(len(taxids1))+" and "+str(len(taxids2))
    all_mrcas=[]
    import operator
    order=operator.itemgetter(-2)
    for taxid1 in taxids1:
        mrcas=[]
        for taxid2 in taxids2:
            mrca,depth1,depth2=find_MRCA(taxid1,taxid2,parents,levels)
            mrcas.append([taxid1,taxid2,mrca,depth1,depth2])
        mrcas=sorted(mrcas,key=order)
        mindepth1=20
        for taxid1,taxid2,mrca,depth1,depth2 in mrcas:
            if depth1<mindepth1:
                mindepth1=depth1
        keep1_mrcas=[]
        for taxid1,taxid2,mrca,depth1,depth2 in mrcas:
            if depth1==mindepth1:
                keep1_mrcas.append([taxid1,taxid2,mrca,depth1,depth2])
        maxdepth2=0
        for taxid1,taxid2,mrca,depth1,depth2 in keep1_mrcas:
            if depth2>maxdepth2:
                maxdepth2=depth2
        keep2_mrcas=[]
        for taxid1,taxid2,mrca,depth1,depth2 in keep1_mrcas:
            if depth2==maxdepth2:
                keep2_mrcas.append([taxid1,taxid2,mrca,depth1,depth2])
        all_mrcas.append(keep2_mrcas)
    return all_mrcas

##### Finds most recent common ancestor - internal function #####
def find_MRCA(taxid1,taxid2,parents,levels):
    reverse=True
    lineage1=get_lineage(taxid1,parents,levels,reverse)
    lineage2=get_lineage(taxid2,parents,levels,reverse)
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

##### Reads in taxid file with non-NCBI synonyms annotated
##### format, name matched to NCBI entry: name|taxid|rank
##### format, name no matched to NCBI entry:
def gather_taxinfo_wsynonyms(taxid_filename,infotypes):
    print "Gathering taxonomic info from file: "+taxid_filename
    taxdict=dict()
    inf=open(taxid_filename)
    lines=inf.read().strip().split('\n')
    inf.close()
    taxannotdict=dict()
    mismatch=dict()
    for i,line in enumerate(lines):
        try:
            taxname,taxid,taxrank=line.strip().split('|')
        except ValueError:
            continue
        if line.startswith('*'):
            taxname=taxname.split('*')[1]
            mismatch_taxname=lines[i-1].strip().split('|')[0]
            try:
                existing=mismatch[taxid]
            except KeyError:
                existing=[]
            existing.append(mismatch_taxname)
            mismatch[taxid]=existing
        taxannotdict[taxid]=[taxname,taxrank]
    print '\t'+str(len(taxannotdict))+' taxa parsed'
    return [taxannotdict,mismatch]

##### Pulls, from the raw iTOL tree newick format file, names of organisms
def pull_org_names_iTOL_newick(itol_infilename):
    inf=open(itol_infilename)
    itoltext=inf.read().strip()
    inf.close()
    import re
    itol_namefinder=re.compile('[(),]{1}[a-zA-Z0-9_.-]+:{1}')
    itol_orgs=itol_namefinder.findall(itoltext)
    itol_orgs=[org[1:-1] for org in itol_orgs]
    itol_orgs=sorted(list(set(itol_orgs)))
    return itol_orgs

##### Creates a mapping from original iTOL org name to simplifed (script) code name
def match_iTOL_to_simplified_names(itol_orgs,taxdict,mismatch):
    itoltaxdict=dict()
    missing=[]
    for itol_org in itol_orgs:
        itol_org_s=simplify_description(itol_org)
        found=False
        for taxid in taxdict.keys():
            taxname=taxdict[taxid][0]
            if itol_org_s==taxname:
                itoltaxdict[taxid]=itol_org
                found=True
        if not found:
            missing.append(itol_org)
    missing2=[]
    for missing_itol_org in missing:
        missing_itol_org_s=simplify_description(missing_itol_org)
        found=False
        for taxid in mismatch.keys():
            taxmismatchnames=mismatch[taxid]
            if missing_itol_org_s in taxmismatchnames:
                itoltaxdict[taxid]=missing_itol_org
                found=True
                break
        if not found:
            missing2.append(missing_itol_org)
    return [itoltaxdict,missing2]

##### Takes the output of the find_MRCA_ext function together with the iTOL names dict generated by
##### match_iTOL_to_simplified_names and outputs a list of iTOL org names that have been mapped to
def output_iTOL_orgs(all_mrcas,itolmm,itoltaxdict,levels):
    print "Selecting final organisms in target taxa"
    itolnames=[]
    for mrcas in all_mrcas:
        found=False
        for mrca in mrcas:
            mrca_taxid=mrca[2]
            try:
                lvl=levels[mrca_taxid]
            except KeyError:
                continue
            if not lvl=='kingdom' and not lvl=='superkingdom' and not lvl=='domain':
                itol_taxid=mrca[1]
                try:
                    itolnamez=itolmm[itol_taxid]
                    itolnames.extend(itolnamez)
                except KeyError:
                    itolname=itoltaxdict[itol_taxid]
                    itolnames.append(itolname)
                found=True
                break
    return itolnames                

##### Simplifies a string (e.g. the iTOL name of an organism from the newick tree file)
##### by removing all non-alphanumeric characters, including underscores, and replacing with
##### whitespace, then rendering lowercase
def simplify_description(description):
    desc=re.sub(r'\W+',' ',description)
    desc=desc.replace('_',' ')
    desc=desc.lower()
    return desc

def simple_getrows(infilename):
    print "Reading in rows from: "+infilename
    inf=open(infilename)
    entries=inf.read().strip().split('\n')
    inf.close()
    print str(len(entries))+" rows parsed"
    return entries

filterset=set([])
parents=make_tax_parents(nodes_filename,filterset)
make_tax_children(node_filename,filterset)
levels=make_tax_level(nodes_filename,filterset)
infotypes=[]
wdir='C:\\Users\\Norman\\Documents\\VCU_research\\Syn3.0'
deg_taxid_filename=wdir+'\\'+'DEG_taxa_wrank.txt'
itol_taxid_filename=wdir+'\\'+'iTOL_taxids_wrank.txt'
degtaxdict,degmm=gather_taxinfo_wsynonyms(deg_taxid_filename,infotypes)
itoltaxdict,itolmm=gather_taxinfo_wsynonyms(itol_taxid_filename,infotypes)
syn3_homologs_taxid_filename=wdir+'\\'+'syn3.0_swissprothomologs_taxids.txt'
syn3_taxids=simple_getrows(syn3_homologs_taxid_filename)
##mrcas=find_MRCA_ext(degtaxdict,itoltaxdict,parents,levels)
allmrcas=find_MRCA_ext(syn3_taxids,itoltaxdict,parents,levels)
itol_newick_infilename=wdir+'\\'+'iTOL_tree_newick.txt'
itol_orgs=pull_org_names_iTOL_newick(itol_newick_infilename)
itoltaxdict2,itolmm2=match_iTOL_to_simplified_names(itol_orgs,itoltaxdict,itolmm)
itol_out_orgs=output_iTOL_orgs(allmrcas,itolmm,itoltaxdict2,levels)
