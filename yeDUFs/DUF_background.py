wdir='C:\\Users\\Norman.Goodacre\\yeDUFs'
neg_kws=['duf','unknown function','uncharacterized protein','uncharacterized domain','uncharacterised protein','uncharacterised domain']

def get_DUFs(DUF_seed_filename):
    inf=open(DUF_seed_filename)
    entries=inf.read().strip().split('#=GF ID')[1:]
    dufinfo=dict()
    for e,entry in enumerate(entries):
        dufac=entry.strip().split('#=GF AC')[1].split('\n')[0].strip()
        try:
            dufid=entry.strip().split('#=GF ID')[1].split('\n')[0].strip()
        except IndexError:
            dufid='NA'
        desc=entry.strip().split('#=GF DE')[1].split('\n')[0].strip()
        isduf=False
        for negkw in neg_kws:
            if negkw in desc.lower():
                isduf=True
        comment_lines=entry.split('#=GF CC')[1:]
        try:
            comment_lastline=comment_lines[-1]
            comment_lastline=comment_lastline.split('\n')[0]
            comment_lines[-1]=comment_lastline
            comment_lines=[line.strip() for line in comment_lines]
            comment=''
            for l,line in enumerate(comment_lines):
                if not l==0:
                    comment+=' '
                comment+=line
            if e==0:
                print comment
            for negkw in neg_kws:
                if negkw in comment.lower():
                    isduf=True
        except IndexError:
            comment='NA'
        if isduf:
            this_dufinfo=[dufac,dufid,desc,comment]
            dufinfo[dufac]=this_dufinfo
    return dufinfo

def get_linkout_proteins(DUF_full_filename,dufs):
    print "Collecting proteins for user-provided set of "+str(len(dufs))+" dufs"
    inf=open(DUF_full_filename)
    pdict=dict()
    c=0
    for i,line in enumerate(inf):
        if i-c==1000000:
            print i
            c=i
        if line.startswith('#=GF AC'):
            pfamid=line.strip().split()[-1].split('.')[0]
        if line.startswith('#=GS'):
            if ' AC ' in line:
                uniprotid=line.strip().split()[-1].strip().split('.')[0]
                try:
                    existing=pdict[pfamid]
                except KeyError:
                    existing=set([])
                existing.add(uniprotid)
                pdict[pfamid]=existing
    return pdict

def yeduf_prots(yedufs,pdict):
    ## store yeduf-containing prots
    yeduf_prots=[]
    for yeduf in list(yedufs):
        ydprots=pdict[yeduf]
        yeduf_prots.extend(list(ydprots))
    yeduf_prots=set(yeduf_prots)
    ## initiate variable for mapping proteins to taxids
    prot_taxdict=dict()
    ## collect yeduf-containing proteins from reference proteomes
    inf=open(wdir+'\\'+'ref_orgs'+'\\'+'reforg_prots.tab')
    yeduf_prots_ref=[]
    for i,line in enumerate(inf):
        sl=line.strip().split('\t')
        prot=sl[0].strip()
        taxid=sl[-1].strip()
        if prot in yeduf_prots:
            yeduf_prots_ref.append(prot)
            prot_taxdict[prot]=taxid
    inf.close()
    ## collect yeduf_containing proteins from non-reference proteomes (x13)
    yeduf_prots_nonref=[]
    import os
    for fn in os.listdir(wdir+'\\'+'nonref_orgs'):
        if not fn.endswith('_prots.csv'):
            continue
        inf=open(wdir+'\\'+'nonref_orgs'+'\\'+fn)
        reader=csv.reader(inf)
        for i,row in enumerate(reader):
            if i==0:
                continue
            prot=row[0].strip()
            taxid=row[-1].strip()
            if prot in yeduf_prots:
                yeduf_prots_nonref.append(prot)
                prot_taxdict[prot]=taxid
    yeduf_prots_filter=[]
    yeduf_prots_filter.extend(yeduf_prots_ref)
    yeduf_prots_filter.extend(yeduf_prots_nonref)
    yeduf_prots_filter=set(yeduf_prots_filter)
    inf.close()
    del reader
    ## write output file, which is in $yeduf $prot $taxid form, per row
    outf=open(wdir+'\\'+'yeDUFs'+'\\'+'yedufprot_yeduf_org.tab','w')
    for yeduf in pdict.keys():
        prots=pdict[yeduf]
        prots_filtered=list(set(prots).intersection(yeduf_prots_filter))
        for prot_f in prots_filtered:
            taxid=prot_taxdict[prot_f]
            outf.write(yeduf+'\t'+prot_f+'\t'+taxid+'\n')
    outf.close()                

def get_taxdict():
    print 'Getting taxdict for yeDUFs, filtering by representative-proteome organisms'
    inf=open(wdir+'\\ref_orgs\\'+'yedufprot_yeduf_reforg.tab')
    taxdict=dict()
    yedufdict=dict()
    for i,line in enumerate(inf):
        prot,yedufs,taxon=line.strip().split('\t')
        yedufs=yedufs.strip().split(',')
        try:
            existing=taxdict[taxon]
        except KeyError:
            existing=[]
        existing.extend(yedufs)
        taxdict[taxon]=existing
        for yeduf in yedufs:
            try:
                existing=yedufdict[yeduf]
            except KeyError:
                existing=set([])
            existing.add(taxon)
            yedufdict[yeduf]=existing
    inf.close()
    return taxdict,yedufdict
                
####### Collects info for all known DUFs
##dufinfo=get_DUFs(DUF_seed_filename)

##### Collects protein accs for user-provided DUF list
##duf_filename='C:\\Users\\Norman.Goodacre\\yeDUFs\\yedufs68.txt'
##duf_inf=open(duf_filename)
##dufs=set(duf_inf.read().strip().split('\n'))
##duf_inf.close()
##DUF_full_filename='C:\\Users\\Norman.Goodacre\\Pfam\\Pfam-A.full'
##pdict=get_linkout_proteins(DUF_full_filename,dufs)

##### Collects yedufs by taxon, where taxon is a reference organism
##taxdict,yedufdict=get_taxdict()
