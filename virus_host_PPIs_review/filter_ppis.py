##########################################################    
##### Begin function and parameter block             #####
##########################################################

### PSI-MS experimental method types that qualify as evidence of "direct AND physical" protein-protein interactions (PPIs)
method_types=['psi-mi:"MI:0397"(two hybrid array)',
'psi-mi:"MI:0096"(pull down)',
'psi-mi:"MI:0018"(two hybrid)',
'psi-mi:"MI:1112"(two hybrid prey pooling approach)',
'psi-mi:"MI:0006"(anti bait coimmunoprecipitation)',
'psi-mi:"MI:0398"(two hybrid pooling approach)',
'psi-mi:"MI:0432"(one hybrid)',
'psi-mi:"MI:0027"(cosedimentation)',
'psi-mi:"MI:1356"(validated two hybrid)',
'psi-mi:"MI:0363"(inferred by author)',
'psi-mi:"MI:2277"(Cr-two hybrid)',
'psi-mi:"MI:1314"(proximity-dependent biotin identification)',
'psi-mi:"MI:0004"(affinity chromatography technology)',
'psi-mi:"MI:0114"(x-ray crystallography)',
'psi-mi:"MI:0071"(molecular sieving)',
'psi-mi:"MI:0030"(cross-linking study)',
'psi-mi:"MI:2222"(inference by socio-affinity scoring)',
'psi-mi:"MI:0424"(protein kinase assay)',
'psi-mi:"MI:0402"(chromatin immunoprecipitation assay)',
'psi-mi:"MI:0111"(dihydrofolate reductase reconstruction)',
'psi-mi:"MI:0029"(cosedimentation through density gradient)',
'psi-mi:"MI:0416"(fluorescence microscopy)',
'psi-mi:"MI:0399"(two hybrid fragment pooling approach)',
'psi-mi:"MI:0019"(coimmunoprecipitation)',
'psi-mi:"MI:0081"(peptide array)',
'psi-mi:"MI:0089"(protein array)',
'psi-mi:"MI:0729"(luminescence based mammalian interactome mapping)',
'psi-mi:"MI:0107"(surface plasmon resonance)',
'psi-mi:"MI:0663"(confocal microscopy)',
'psi-mi:"MI:0053"(fluorescence polarization spectroscopy)',
'psi-mi:"MI:0809"(bimolecular fluorescence complementation)',
'psi-mi:"MI:0400"(affinity technology)',
'psi-mi:"MI:0065"(isothermal titration calorimetry)',
'psi-mi:"MI:0112"(ubiquitin reconstruction)',
'psi-mi:"MI:0028"(cosedimentation in solution)',
'psi-mi:"MI:0813"(proximity ligation assay)',
'psi-mi:"MI:0276"(blue native page)',
'psi-mi:"MI:0077"(nuclear magnetic resonance)']
method_types_short=[m.split('(')[1].split(')')[0] for m in method_types]
method_types_code=[m.split('"')[1].split('"')[0] for m in method_types]
method_types_biogrid=['Two-hybrid','Affinity Capture-Western','Far Western','Co-crystal Structure',
                      'Protein-peptide','Affinity Capture-Luminescence','Proximity Label-MS',
                      'FRET','Co-purification','Affinity Capture-MS']
### PSI-MS interaction types that qualify as evidence of "direct AND physical" protein-protein interactions (PPIs)
### ..so far, only one
int_types=['psi-mi:"MI:0407"(direct interaction)']
int_types_short=[i.split('(')[1].split(')')[0] for i in int_types]
int_types_code=[i.split('"')[1].split('"')[0] for i in method_types]
import re
mifinder=re.compile('MI:[0-9]{4}')
ufinder=re.compile('[A-Z]{1}[A-Z0-9]{5}[A-Z0-9]*')

##### BioGRID lists PPIs by their gene ids; this requires mapping to uniprot ids. List gene id -> uniprot .tab mapping file(s)
def get_idmap(ppis_dir,idmap_filenames):
    idmap=dict()
    for idmap_filename in idmap_filenames:
        inf=open(ppis_dir+'\\'+idmap_filename)
        for line in inf:
            id1,uniprotid=line.strip().split('\t')
            try:
                existing=idmap[id1]
            except KeyError:
                existing=[]
            existing.append(uniprotid)
            idmap[id1]=existing
        inf.close()
    return idmap

##### Fetches direct, physical PPIs from a file containing PPIs in PSI-MS format
##### the method types and interaction type(s) are listed above ^
##### Note that for some PSI-MS databases, uniprot IDs aren't used, so a mapping parameter
##### is optional (idmap)
def fetch_directinteractions_psims(ppis_dir,ppis_filename,method_types,int_types,idmap):
    print "Parsing DIRECT PPIs from: "+ppis_filename
    out_ppis=[]
    inf=open(ppis_dir+'\\'+ppis_filename)
    outf=open(ppis_dir+'\\direct_PPIS\\'+ppis_filename+'.DIRECT.tab','w')
    c=0 # tracks the number of PPIs total that are parsed
    d=0 # tracks the number of PPIs that pass the exp method and int type filter
    e=0 # tracks the number of PPIs that pass the exp method and int type filter, as well as have two uniprot IDs (present or by mapping)
    for i,line in enumerate(inf):
        if i==0:
            continue
        sl=line.split('\t')
        imethods=mifinder.findall(sl[6])
        itypes=mifinder.findall(sl[11])
        c+=1
        valid=False
        for imethod in imethods:
            if imethod in method_types_code:
                valid=True
        for itype in itypes:
            if itype in int_types_code:
                valid=True
        if valid:
            d+=1
        else:
            continue
        ppis=parse_psims_foruniprotids(line,idmap,False)
        if ppis=='none':
            continue
        e+=1
        out_ppis.extend(ppis)
        for ppi in ppis:
            outf.write('\t'.join(ppi)+'\n')
    print str(e)+' interactions with uniprotids out of '+str(d)+' that are experimentally valid out of '+str(c)+' total parsed'
    return out_ppis
    outf.close()

##### Fetches direct, physical PPIs from BioGRID
##### the method types and interaction type(s) are listed above ^ but are the "short" version (name only, without identifiers)
##### Note that BioGRID uses geneids and BioGRID ids, so a mapping file is mandatory (-> uniprot ids)
def fetch_directinteractions_biogrid(ppis_dir,ppis_filename,method_types_biogrid,idmap):
    print "Parsing DIRECT PPIs from: "+ppis_filename
    out_ppis=[]
    inf=open(ppis_dir+'\\'+ppis_filename)
    outf=open(ppis_dir+'\\direct_PPIS\\'+ppis_filename+'.DIRECT.tab','w')
    c=0 # tracks the number of PPIs total that are parsed
    d=0 # tracks the number of PPIs that pass the exp method and int type filter
    e=0 # tracks the number of PPIs that pass the exp method and int type filter, as well as have two uniprot IDs (present or by mapping)
    for i,line in enumerate(inf):
        if i==0:
            continue
        sl=line.split('\t')
        imethod=sl[11]
        itype=sl[12]
        c+=1
        if imethod in method_types_biogrid:
            valid=True
            d+=1
        else:
            valid=False
            continue
        ppis=parse_psims_foruniprotids(line,idmap,True)
        if ppis=='none':
            continue
        e+=1
        out_ppis.extend(ppis)
        for ppi in ppis:
            outf.write('\t'.join(ppi)+'\n')
    print str(e)+' interactions with uniprotids out of '+str(d)+' that are experimentally valid out of '+str(c)+' total parsed'
    return out_ppis
    outf.close()

def parse_psims_foruniprotids(line,idmap,biogrid):
    sl=line.strip().split('\t')
    if biogrid:
        a1,b1,a2,b2=sl[1:5]
    else:
        a1,b1,a2,b2=sl[:4]
    a=a1+'|'+a2
    b=b1+'|'+b2
    uniprota=[]
    uniprotb=[]
    uniprota.extend(ufinder.findall(a))
    uniprotb.extend(ufinder.findall(b))
    if len(uniprota)>0 and len(uniprotb)>0:
        ints=[]
        for ua in uniprota:
            for ub in uniprotb:
                ints.append([ua,ub])
        return ints
    else:
        ints=parse_psims_foruniprotids2(a,b,idmap,biogrid)
        return ints

def parse_psims_foruniprotids2(a,b,idmap,biogrid):
    uniprota2=[]
    uniprotb2=[]
    if biogrid:
        ida=a.strip().split('|')
        idb=b.strip().split('|')
    else:
        ida=[]
        for i in a.strip().split('|'):
            try:
                ida.append(i.split(':')[1])
            except IndexError:
                continue
        idb=[]
        for i in b.strip().split('|'):
            try:
                idb.append(i.split(':')[1])
            except IndexError:
                continue
    for i in ida:
        try:
            u=idmap[i]
            uniprota2.extend(u)
        except KeyError:
            continue
    for i in idb:
        try:
            u=idmap[i]
            uniprotb2.extend(u)
        except KeyError:
            continue
    if len(uniprota2)>0 and len(uniprotb2)>0:
        ints=[]
        for ua in uniprota2:
            for ub in uniprotb2:
                ints.append([ua,ub])
        return ints
    else:
        return 'none'

##########################################################    
##### End function and parameter block               #####
##########################################################



##########################################################    
##### Begin execute block                            #####
##########################################################
ppis_dir='E:\\PPI_dbs'
intact_filename='intact_03282018.txt'
dip_filename='dip20170205.txt'
mint_filename='MINT_03282018_taxid9606.txt' # the MINT file is for humans only, since MINT has an inefficient page-display form rather than download link
biogrid_filename='BIOGRID-ALL-3.4.158.tab2.txt'
matrixdb_filename='matrixdb_171107_FULL27.mitab'
idmap_filename='biogrid_geneids_03282018_mapping.tab'
idmap=get_idmap(ppis_dir,[idmap_filename])
##intact_ints=fetch_directinteractions_psims(ppis_dir,intact_filename,method_types,int_types,idmap)
##dip_ints=fetch_directinteractions_psims(ppis_dir,dip_filename,method_types,int_types,idmap)
##mint_ints=fetch_directinteractions_psims(ppis_dir,mint_filename,method_types,int_types,idmap)
##matrix_ints=fetch_directinteractions_psims(ppis_dir,matrixdb_filename,method_types,int_types,idmap)
biogrid_ints=fetch_directinteractions_biogrid(ppis_dir,biogrid_filename,method_types_biogrid,idmap)
##########################################################    
##### End execute block                              #####
##########################################################
