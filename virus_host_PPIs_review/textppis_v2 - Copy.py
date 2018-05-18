import re
import csv

def get_synonyms_human(uniprot_dl_filename):
    inf=open(uniprot_dl_filename)
    reader=csv.reader(inf)
    synonyms=dict()
    for i,row in enumerate(reader):
        uniprotac=row[0]
        names=row[8]
        names=names.strip().split('(')
        names=[name.split(')')[0].strip() for name in names]
        names=[simplify_string(name) for name in names]
        keepnames=[]
        for n,name in enumerate(names):
            name=name.strip()
            if name in neg_human_kws:
                print name
                continue
            if len(name)>=3:
                primary_name=name
                break
        for name in names[n:]:
            name=name.strip()
            if name in neg_human_kws:
                continue
            if len(name)<3:
                continue
            keepnames.append(name)
        keepnames=[' '+keepname+' ' for keepname in keepnames]
        synonyms[primary_name]=keepnames
    inf.close()
    del reader
    return synonyms

def get_synonyms_viral(viral_namefile):
    synonyms=dict()
    vnamesall=[]
    inf=open(viral_namefile)
    reader=csv.reader(inf)
    for i,row in enumerate(reader):
        if i==0:
            continue
        vnames=row[0].strip().split(',')
        vnames_keep=[]
        for vname in vnames:
            vname=vname.strip()
            if len(vname)>=3:
                vnames_keep.append(' '+vname+' ')
        if i==0:
            continue
        vnamesall.extend(vnames_keep)
        if len(row[1])>0:
            primary_name=' '+row[1].strip()+' '
            for vname in vnames_keep:
                synonyms[vname]=primary_name
    del reader
    inf.close()
    vnamesall=list(set(vnamesall))
    return [vnamesall,synonyms]

def scan_abstract_ctrl(abstracts_filename,synonyms_human,names_hiv,synonyms_hiv,out_logfilename):
    inf=open(abstracts_filename)
    abstracts=inf.read().strip().split('PMID: ')
    inf.close()
    abs_ppis_main=[]
    c=0
    out_logfile=open(out_logfilename,'w')
    ppis_abstract_freq=dict()
    print str(len(abstracts))+' abstracts in total'
    for a,abstract in enumerate(abstracts):
        if a-c==100:
            print str(a)+' abstracts scanned for virus-human PPIs'
            c=a
        abstract_s=simplify_string(abstract)
        abs_ppis,intexts,ppi_abstract_freq=scan_abstract(abstract_s,synonyms_human,names_hiv,synonyms_hiv,ppis_abstract_freq)
        abs_ppis_main.extend(abs_ppis)
        out_logfile.write('\n'.join(intexts)+'\n\n')
    out_logfile.close()
    return abs_ppis_main,ppis_abstract_freq

def scan_abstract(abstract_s,synonyms_human,names_hiv,synonyms_hiv,ppis_abstract_freq):
    abs_ppis=[]
    hiv_intexts=[]
    intexts=[]
    for name_hiv in names_hiv:
        if name_hiv in abstract_s:
            hivpos=abstract_s.find(name_hiv)
            hiv_intext=abstract_s[hivpos-25:hivpos+30]
            try:
                primary_name_hiv=synonyms_hiv[name_hiv]
            except KeyError:
                primary_name_hiv=name_hiv
            for primary_name_human in synonyms_human.keys():
                names_human=synonyms_human[primary_name_human]
                for name_human in names_human:
                    if name_human in abstract_s:
                        #####
                        humanpos=abstract_s.find(name_human)
                        human_intext=abstract_s[humanpos-25:humanpos+30]
                        intexts.append(name_hiv+'...'+hiv_intext+'... ...'+name_human+'...'+human_intext)
                        abs_ppis.append([primary_name_hiv,primary_name_human])
    found=set([])
    for abs_ppi in abs_ppis:
        tag='__'.join(abs_ppi)
        if tag in found:
            continue
        try:
            existing=ppis_abstract_freq[tag]
        except KeyError:
            existing=0
        existing+=1
        ppis_abstract_freq[tag]=existing
        found.add(tag)
    return [abs_ppis,intexts,ppis_abstract_freq]

def abs_ppis_nr(abs_ppis):
    found=set([])
    abs_ppis_nr=[]
    for abs_ppi in abs_ppis:
        if '_'.join(abs_ppi) in found:
            continue
        else:
            abs_ppis_nr.append(abs_ppi)
            found.add('_'.join(abs_ppi))
    return abs_ppis_nr

neg_human_kws=['cytosine','dna','trna','fragment','alpha','beta','protein can','oct 1','1 6','poly','van','aid','set','delta']

import sys
sys.path.append('C:\\Users\\Norman.Goodacre\\Documents\\TOOLBOX')
from text_functions import *
uniprot_dl_filename='C:\\Users\\Norman.Goodacre\\Documents\\virus_ppis\\human_prots_info_wHIV1.csv'
viral_namefile='C:\\Users\\Norman.Goodacre\\Documents\\virus_ppis\\hiv1_prot_names3.csv'
abstracts_filename='C:\\Users\\Norman.Goodacre\\Documents\\virus_ppis\\pubmed_HIV-human PPIs abstracts.txt'
synonyms_human=get_synonyms_human(uniprot_dl_filename)
names_hiv,synonyms_hiv=get_synonyms_viral(viral_namefile)
out_logfilename='C:\\Users\\Norman.Goodacre\\Documents\\virus_ppis\\textppis_spotcheck.txt'
ppis_abstract,ppis_abstract_freq=scan_abstract_ctrl(abstracts_filename,synonyms_human,names_hiv,synonyms_hiv,out_logfilename)
ppis_abstract_nr=abs_ppis_nr(ppis_abstract)
