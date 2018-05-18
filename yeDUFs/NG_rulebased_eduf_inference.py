main_dir='C:\\Users\\Norman\\Documents\\VCU_research'
DEG_main_dir=main_dir+'\\'+'DEG10'
euk_essential_folder='deg-e-10'
euk_nonessential_folder='deg-ne-10'
euk_ne_path=DEG_main_dir+'\\'+euk_nonessential_folder
euk_e_path=DEG_main_dir+'\\'+euk_essential_folder
DUFsfn='DUFs_Pfamv27.txt'
eduf_proteins_flatfn='yeCase1Case2DUFs_proteins.txt'
cerevis_e_pfamfn='559292_e_pfams.txt'
pombe_e_pfamfn='284812_e_pfams.txt'
pombe_ne_pfamfn='284812_ne_pfams.txt'

def get_DUFs():
    infile=open(main_dir+'\\'+DUFsfn)
    dufs=infile.read().strip().split('\n')
    dufs=set(dufs)
    return dufs
    
def case1(e_pfams_fn):
    case1_hits=[]
    inf=open(e_pfams_fn)
    for line in inf:
        try:
            prot,pfams=line.strip().split()
        except ValueError:
            continue
        pfams=set(pfams.strip().split(','))
        if len(pfams.intersection(dufs))==len(pfams):
        ##if len(pfams.difference(dufs))==0:
            out=[prot]
            out.extend(list(pfams))
            case1_hits.append(out)
    return case1_hits

def case2(e_pfams_fn,ne_pfams_fn):
    e_arch=case2_dictmaker(e_pfams_fn)
    ne_arch=case2_dictmaker(ne_pfams_fn)
    print len(e_arch)
    print len(ne_arch)
    case2_hits=[]
    for eprot in e_arch.keys():
        hit=False
        e_arch_=e_arch[eprot]
        e_arch_dufs=set(e_arch_).intersection(dufs)
        if len(e_arch_dufs)==0:
            continue
        for neprot in ne_arch.keys():
            ne_arch_=set(ne_arch[neprot])
            ne_arch_diff=ne_arch_.difference(set(e_arch_))
            if ne_arch_diff==e_arch_dufs:
                hit=True
##                case2_hit=[eprot,neprot]
##                case2_hit.extend(e_arch_dufs)
##                case2_hits.append(case2_hit)
        if hit:
            case2_hit=[eprot]
            case2_hit.extend(list(e_arch_dufs))
            case2_hits.append(case2_hit)
    return case2_hits
            

def case2_dictmaker(pfams_fn):
    inf=open(pfams_fn)
    arch=dict()
    for line in inf:
        try:
            prot,pfams=line.strip().split()
        except ValueError:
            continue
        pfams=pfams.strip().split(',')
        if len(pfams)==0:
            continue
        arch[prot]=pfams
    return arch

def case3(e_pfams_fn,ne_pfams_fn):
    e_arch=case2_dictmaker(e_pfams_fn)
    ne_arch=case2_dictmaker(ne_pfams_fn)
    e_rep=case3_inverse_dictmaker(e_arch)
    ne_rep=case3_inverse_dictmaker(ne_arch)
    case3_hits=[]
    for epfam in e_rep.keys():
        e_prots=e_rep[epfam]
        if len(e_prots)<=1:
            continue
        if not epfam in dufs:
            continue
        if not epfam in ne_rep.keys():
            case3_hit=[epfam]
            case3_hit.extend(e_prots)
            case3_hits.append(case3_hit)
    return case3_hits

def case3_inverse_dictmaker(arch):
    pfam_rep=dict()
    for prot in arch.keys():
        pfams=arch[prot]
        for pfam in pfams:
            try:
                existing=pfam_rep[pfam]
            except KeyError:
                existing=[]
            existing.append(prot)
            pfam_rep[pfam]=existing
    return pfam_rep

def get_description_proteinflatfile(flatfn,allprots):
    desc=dict()
    inf=open(main_dir+'\\'+flatfn)
    entries=inf.read().strip().split('\n//\n')
    inf.close()
    for entry in entries:
        alldesc=entry.strip().split('RecName:')[1].split('AltName:')[0].strip()
        fulldesc=alldesc.split('Full=')[1].strip().split(';')[0]
        lines=entry.strip().split('\n')
        for line in lines:
            if line.startswith('AC '):
                acs=line.strip().split()[1:]
                acs=[ac.strip().split(';')[0] for ac in acs]
                acs=set(acs)
                try:
                    ac=list(acs.intersection(allprots))[0]
                except IndexError:
                    print acs
        desc[ac]=fulldesc
    return desc

def get_eDUF(dufs,pfamsfn):
    edufs=dict()
    inf=open(euk_e_path+'\\'+pfamsfn)
    for line in inf:
        try:
            prot,pfams=line.strip().split()
        except ValueError:
            continue
        pfams=pfams.strip().split(',')
        edufs_=set(pfams).intersection(dufs)
        pfams_=set(pfams).difference(dufs)
        edufs[prot]=[pfams_,edufs_]
    return edufs
   

dufs=get_DUFs()
case1_cerevis=case1(euk_e_path+'\\'+cerevis_e_pfamfn)
case1_pombe=case1(euk_e_path+'\\'+pombe_e_pfamfn)
case2_cerevis=case2(euk_e_path+'\\'+cerevis_e_pfamfn,euk_ne_path+'\\'+pombe_ne_pfamfn)
case2_pombe=case2(euk_e_path+'\\'+pombe_e_pfamfn,euk_ne_path+'\\'+pombe_ne_pfamfn)
allprots=[]
allprots.extend([c[0] for c in case1_cerevis])
allprots.extend([c[0] for c in case1_pombe])
allprots.extend([c[0] for c in case2_cerevis])
allprots.extend([c[0] for c in case2_pombe])
allprots=set(allprots)
desc=get_description_proteinflatfile(eduf_proteins_flatfn,allprots)
print len(dufs)
case3_cerevis=case3(euk_e_path+'\\'+cerevis_e_pfamfn,euk_ne_path+'\\'+pombe_ne_pfamfn)
case3_pombe=case3(euk_e_path+'\\'+pombe_e_pfamfn,euk_ne_path+'\\'+pombe_ne_pfamfn)
pfaminfo_cerevis=get_eDUF(dufs,cerevis_e_pfamfn)
pfaminfo_pombe=get_eDUF(dufs,pombe_e_pfamfn)
print 'Protein full names, case 1 S. cerevisiae: '
for protinf in case1_cerevis:
    ac=protinf[0]
    pfams,dufs=pfaminfo_cerevis[ac]
    print ac+';'+','.join(pfams)+';'+','.join(dufs)+';'+'_'.join(desc[ac].strip().split())
print '\n'
print 'Protein full names, case 1 S. pombe: '
for protinf in case1_pombe:
    ac=protinf[0]
    pfams,dufs=pfaminfo_pombe[ac]
    print ac+';'+','.join(pfams)+';'+','.join(dufs)+';'+'_'.join(desc[ac].strip().split())
print '\n'
print 'Protein full names, case 2 S. cerevisiae: '
for protinf in case2_cerevis:
    ac=protinf[0]
    pfams,dufs=pfaminfo_cerevis[ac]
    print ac+';'+','.join(pfams)+';'+','.join(dufs)+';'+'_'.join(desc[ac].strip().split())
print '\n'
print 'Protein full names, case 2 S. pombe: '
for protinf in case2_pombe:
    ac=protinf[0]
    pfams,dufs=pfaminfo_pombe[ac]
    print ac+';'+','.join(pfams)+';'+','.join(dufs)+';'+'_'.join(desc[ac].strip().split())
print '\n'
