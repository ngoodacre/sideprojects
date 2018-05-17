taxdir='F:\\TAXONOMY\\April2017'
gi_map_filename=taxdir+'\\'+'gi_taxid_nucl.dmp'
names_filename=taxdir+'\\'+'taxdump.tar'+'\\'+'names.dmp'
##nodes_filename=taxdir+'\\'+'taxdump.tar'+'\\'+'nodes.dmp'
ICTV_master_filename='F:\\TAXONOMY'+'\\'+'ICTV_master_sheet.csv'
ICTV_orgs_filename='F:\\TAXONOMY'+'\\'+'ICTV_viral_species.txt'
ICTV_NCBI_par_filename='F:\\TAXONOMY'+'\\'+'ICTV_NCBI_viral_taxa.csv'

##### Removes any PPIs that have not been identified as direct, physical interactions using a set of pre-extracted filters from 
def filter_ppis_directbinding(

def collect_virusmentha_ppis(virusmentha_ppifile,reforg,targetlevel):
    print 'Collecting virusmentha PPIs, grouped by virus family'
    inf=open(virusmentha_ppifile)
    ppis_byorg=dict()
    for i,line in enumerate(inf):
        if i==0:
            continue
        sl=line.strip().split(';')
        try:
            prot1,gene1,tax1,fam1,prot2,gene2,tax2,fam2,score,pmid=sl
        except ValueError:
            print i
            continue
        if tax1==reforg:
            outprot1=prot1
            outgene1=gene1
            outtax1=tax1
            outfam1=fam1
            outprot2=prot2
            outgene2=gene2
            outtax2=tax2
            outfam2=fam2
        else:
            outprot1=prot2
            outgene1=gene2
            outtax1=tax2
            outfam1=fam2
            outprot2=prot1
            outgene2=gene1
            outtax2=tax1
            outfam2=fam1
        try:
            if targetlevel=='family':
                existing=ppis_byorg[outfam2]
            if targetlevel=='species':
                existing=ppis_byorg[outtax2]
        except KeyError:
            existing=[]
        existing.append([outprot1,outprot2])
        if targetlevel=='family':
            ppis_byorg[outfam2]=existing
        if targetlevel=='species':
            ppis_byorg[tax2]=existing
    return ppis_byorg

def collect_taxinfo(taxinfo_filename,filterset):
    print "Collecting taxonomic information from file: "+taxinfo_filename
    import re
    taxidfinder=re.compile('NCBI_TaxID=[0-9]+[;\w]{1}')
    inf=open(taxinfo_filename)
    taxdict=dict()
    for i,line in enumerate(inf):
        if line.startswith('AC   '):
            uniprotacs=line.strip().split()[1:]
            uniprotacs=[uniprotac.strip().split(';')[0] for uniprotac in uniprotacs]
        if line.startswith('OX   '):
            try:
                taxid=taxidfinder.findall(line.strip())[0]
            except IndexError:
                print line
                print taxidfinder.findall(line.strip())
            taxid=taxid.split('=')[1][:-1]
            for uniprotac in uniprotacs:
                if uniprotac in filterset:
                    taxdict[uniprotac]=taxid
    inf.close()
    return taxdict

def collect_virhostnet_hpidb_ppis(db_filename,viral_taxdict,humanproteins,sp_to_fam,map_up,stage):
    if stage==1:
        print 'Collecting viral and human proteins from '+db_filename
    if stage==2:
        print 'Collecting viral and human PPIs, grouped by virus, from '+db_filename
    inf=open(db_filename)
    if stage==1:
        filterprots=set([])
    if stage==2:
        ppis_byorg=dict()
    for i,line in enumerate(inf):
        sl=line.strip().split()
        if 'virhostnet' in db_filename.lower():
            uniprotac1=sl[0].split(':')[1].strip()
            uniprotac2=sl[1].split(':')[1].strip()
        if 'hpidb' in db_filename.lower():
            if i==0:
                continue
            uniprotac1=sl[0].split(':')[-1].strip()
            uniprotac2=sl[3].split(':')[-1].strip()
        viral1=False
        viral2=False
        human1=False
        human2=False
        if stage==1:
            filterprots.add(uniprotac1)
            filterprots.add(uniprotac2)
        if stage==2:
            try:
                taxid1=viral_taxdict[uniprotac1]
                if mapup:
                    taxid1=sp_to_fam[taxid1]
                viral1=True
            except KeyError:
                if uniprotac1 in humanproteins:
                    human1=True
            try:
                taxid2=viral_taxdict[uniprotac2]
                if mapup:
                    taxid2=sp_to_fam[taxid2]
                viral2=True
            except KeyError:
                if uniprotac2 in humanproteins:
                    human2=True
            if viral1 and not viral2:
                try:
                    existing=ppis_byorg[taxid1]
                except KeyError:
                    existing=[]
                existing.append([uniprotac2,uniprotac1])
                ppis_byorg[taxid1]=existing
            if viral2 and not viral1:
                try:
                    existing=ppis_byorg[taxid2]
                except KeyError:
                    existing=[]
                existing.append([uniprotac1,uniprotac2])
                ppis_byorg[taxid2]=existing
    if stage==1:
        return filterprots
    if stage==2:
        return ppis_byorg


                        
reforg='9606'
virusmentha_ppifilename='E:\\virus_host_ppis\\virusmentha_9606'
virhostnet_ppifilename='E:\\virus_host_ppis\\virhostnet2.0.txt'
hpidb_ppifilename='E:\\virus_host_ppis\\hpidb2.tab'
import sys
sys.path.append('E:\\TOOLBOX')
from basic_tax_functions import *
from sequence_record_functions import combine_dicts
from sequence_record_functions import combine_sets
from sequence_record_functions import make_dict
from sequence_record_functions import get_gis_flatfile
from sequence_record_functions import combine_dict_listoflists

tax2name=map_ncbi_taxids_to_taxnames(set([]),False,False)
viral_sprot_filename='E:\\virus_host_ppis\\uniprot_sprot_viruses.dat'
viral_trembl_filename='E:\\virus_host_ppis\\uniprot_trembl_viruses.dat'
human_sprot_filename='E:\\virus_host_ppis\\uniprot_sprot_human.dat'
human_trembl_filename='E:\\virus_host_ppis\\uniprot_trembl_human.dat'

##stage=1
##filterprots_virhostnet=collect_virhostnet_hpidb_ppis(virhostnet_ppifilename,dict(),dict(),dict,False,stage)
##filterprots_hpidb=collect_virhostnet_hpidb_ppis(hpidb_ppifilename,dict(),dict(),dict,False,stage)
##filterprots=combine_sets(filterprots_virhostnet,filterprots_hpidb)
##viral_taxdict_sprot=collect_taxinfo(viral_sprot_filename,filterprots)
##viral_taxdict_trembl=collect_taxinfo(viral_trembl_filename,filterprots)
##human_taxdict_sprot=collect_taxinfo(human_sprot_filename,filterprots)
##human_taxdict_trembl=collect_taxinfo(human_trembl_filename,filterprots)
##islist=False
##viral_taxdict=combine_dicts(viral_taxdict_sprot,viral_taxdict_trembl,islist)
##human_taxdict=combine_dicts(human_taxdict_sprot,human_taxdict_trembl,islist)

viral_mapping_filename='E:\\virus_host_ppis\\virus_taxonomic_mapping.txt'
viral_taxdict=make_dict(viral_mapping_filename,'\t',False)
humanproteins_filename='E:\\virus_host_ppis\\humanproteins_virushostnet_hpidb.txt'
humanproteins=get_gis_flatfile(humanproteins_filename)
parents=make_tax_parents(nodes_filename,dict())
levels=make_tax_level(nodes_filename,dict())
target_level='family'
asdict=True
viral_sp_to_fam=climb_tax_tree_main(viral_taxdict.values(),parents,levels,target_level,asdict)
target_level='species'
viral_sp_to_sp=climb_tax_tree_main(viral_taxdict.values(),parents,levels,target_level,asdict)

stage=2
mapup=True
ppis_byorg_virusmentha=collect_virusmentha_ppis(virusmentha_ppifilename,reforg,'family')
ppis_byorg_virhostnet=collect_virhostnet_hpidb_ppis(virhostnet_ppifilename,viral_taxdict,humanproteins,viral_sp_to_fam,mapup,stage)
ppis_byorg_hpidb=collect_virhostnet_hpidb_ppis(hpidb_ppifilename,viral_taxdict,humanproteins,viral_sp_to_fam,mapup,stage)
ppis_byorg_virhostnet_hpidb=combine_dict_listoflists(ppis_byorg_virhostnet,ppis_byorg_hpidb)
ppis_byorg=combine_dict_listoflists(ppis_byorg_virhostnet_hpidb,ppis_byorg_virusmentha)

stage=2
mapup=True
ppis_byorg_virusmentha=collect_virusmentha_ppis(virusmentha_ppifilename,reforg,'species')
ppis_byorg_virhostnet_sp=collect_virhostnet_hpidb_ppis(virhostnet_ppifilename,viral_taxdict,humanproteins,viral_sp_to_sp,mapup,stage)
ppis_byorg_hpidb_sp=collect_virhostnet_hpidb_ppis(hpidb_ppifilename,viral_taxdict,humanproteins,viral_sp_to_sp,mapup,stage)
ppis_byorg_virhostnet_hpidb_sp=combine_dict_listoflists(ppis_byorg_virhostnet_sp,ppis_byorg_hpidb_sp)
ppis_byorg_sp=combine_dict_listoflists(ppis_byorg_virhostnet_hpidb_sp,ppis_byorg_virusmentha)
##for taxid in ppis_byorg.keys():
##    try:
##        print tax2name[taxid]+'\t'+str(len(ppis_byorg[taxid]))
##    except KeyError:
##        print taxid+' NOT FOUND'
    
