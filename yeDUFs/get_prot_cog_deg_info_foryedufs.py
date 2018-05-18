

import csv


def get_prots_pfams(org_prot_filename,reverse):
    inf=open(org_prot_filename)
    reader=csv.reader(inf)
    d=dict()
    for i,row in enumerate(reader):
        prot=row[0].strip()
        pfaminfo=row[4].strip().split(';')
        if len(pfaminfo)>0:
            pfams=pfaminfo[:-1]
            if reverse:
                for pfam in pfams:
                    pfam=pfam.strip()
                    try:
                        existing=d[pfam]
                    except KeyError:
                        existing=[]
                    existing.append(prot)
                    d[pfam]=existing
            else:
                d[prot]=pfams
        else:
            continue
    return d

def get_cogs_prots(cog_filename,protfilter,reverse):
    d=dict()
    inf=open(cog_filename)
    for i,line in enumerate(inf):
        sl=line.strip().split('\t')
        cogid=sl[1].strip()
        protinfosfull=sl[5].strip()
        protinfos=[protinfofull.strip().split('.')[:2] for protinfofull in protinfosfull]
        if reverse:
            for protinfo in protinfos:
                try:
                    existing=d[protinfo[1]]
                except KeyError:
                    existing=[]
                existing.append([cogid,protinfo[0]])
                d[protinfo[1]]=existing
        else:
            d[cogid]=protinfos
        return d

wdir='F:\\yeDUFs\\April2017'
org_prot_filename=wdir+'\\'+'Sc_uniprot_info_042017.csv'
d=get_prots_pfams(org_prot_filename,True)
