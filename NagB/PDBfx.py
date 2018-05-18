import os
import sys
docking_dir='F:\\SNPs@interface\\Docking'
#os.chdir('F:\\SNPs@interface\\Docking\\HADDOCK')
#sys.path.append('C:\\BCHB524\\Python25\\Lib\\site-packages')
import Bio.PDB.PDBParser
parser=Bio.PDB.PDBParser()

# imports pdb structures. All you need to pass as a parameter is a list of
# pdb structure names
def import_pdb_structures(pdb_filenames_list):
    pdb_structures=dict()
    for i,pdb_filename in enumerate(pdb_filenames_list):
        # This is where the pdb structures get their names (it's the filename sans '.pdb' extension)
        pdb_name=pdb_filename.split('.')[0]
        structure=parser.get_structure(pdb_name,pdb_filename)
        pdb_structures[pdb_filename.split('.')[0]]=structure
    return pdb_structures

# Will take your imported structures (param 1) and extract the residues present
# in the chains you pass as param 2. Param 2 is a dictionary where the keys
# should correspond to the keys in pdb_structures (i.e. the name of the pdb files
# but without the '.pdb' extension) and the values are lists specifying
# which chains you're interested in
# Ex. param 2: {'1KFX':['L','S'], '1KFU':['L','S']}
def make_residues_dict(pdb_structures,pdb_chains):
    pdb_resdict=dict()
    for pdb_name in pdb_structures.keys():
        structure=pdb_structures[pdb_name]
        model=structure[0]
        resdict=dict()
        chainids=pdb_chains[pdb_name]
        for chainid in chainids:
            chain=model[chainid]
            for res in chain:
                resdict[chainid+'_'+str(res.get_id()[1])]=res
        pdb_resdict[pdb_name]=resdict
    return pdb_resdict

# We will be aligning by alpha-carbons below. Also, we need to check that
# the two atom lists to be aligned are of the identical length (all residue positions
# must be shared)
def get_fixed_moving_CAlists(resdict1,resdict2,asdict,chainmap,offset):
    fixed=[]
    moving=[]
    if asdict:
        fixed=dict()
        moving=dict()
    for key in resdict1.keys():
        try:
            chain1,pos1=key.strip().split('_')
            chain2=chainmap[chain1]
            pos2=str(int(pos1)+offset)
            key2=chain2+'_'+pos2
            res2=resdict2[key2]
        except KeyError:
            continue
        res1=resdict1[key]
        if not asdict:
            try:
                fixed.append(res1['CA'])
                moving.append(res2['CA'])
            except KeyError:
                continue
        else:
            try:
                fixed[key]=res1['CA']
                moving[key]=res2['CA']
            except KeyError:
                continue
    return [fixed,moving]
    
# Unfortunately rmsd calculation can only be done pairs at a time (as far as I know),
# so here, submit a pair of residue dictionaries, each corresponding to a
# pdb file, at a time.
def get_rms(fixed,moving):
    superimposer=Bio.PDB.Superimposer()
    superimposer.set_atoms(fixed,moving)
    return superimposer.rms

# Takes all your input pdbs and return the name of the one with the lowest cumulative rmsd with all
# the others. This it does in 'best_pdb' mode. In 'rmsd' mode, will return the dictionary of dictionaries
# (rmsd_dict) containing all rmsd values
def find_representative_pdb(pdb_resdict,mode):
    rmsd_dict=dict()
    best_rmsd=float(1000)
    for pdb_name1 in pdb_resdict.keys():
        rmsd_sum=float(0)
        rmsd_dict_sub=dict()
        for pdb_name2 in pdb_resdict.keys():
            if pdb_name2==pdb_name1:
                continue
            fixed,moving=get_fixed_moving_CAlists(pdb_resdict[pdb_name1],pdb_resdict[pdb_name2],False)
            rms=get_rms(fixed,moving)
            rmsd_dict_sub[pdb_name2]=rms
            rmsd_sum+=rms
        rmsd_dict[pdb_name1]=rmsd_dict_sub
        print pdb_name1+' '+str(rmsd_sum)
        if rmsd_sum < best_rmsd:
            best_pdb=pdb_name1
            best_rmsd=rmsd_sum
    if mode=='best_pdb':
        return best_pdb
    if mode=='rmsd':
        return rmsd_dict

# fixed_name : pdb structure (or it's key in the pdb_resdict dictionary, rather)
# you want to align to. All others align to this one.
# The method find_representative_pdb, below, will allow you to choose the pdb structure that aligns best,
# on average, with all of the other structures. You wouldn't want to align to some
# shitty structure. 
def align_all_to_common_structure(fixed_name,pdb_resdict,pdb_structures,fixed_chain,pdb_chains):
    io=Bio.PDB.PDBIO()
    superimposer=Bio.PDB.Superimposer()
    for pdb_name in pdb_resdict.keys():
        if pdb_name==fixed_name:
            continue
        fixed,moving=get_fixed_moving_CAlists(pdb_resdict[fixed_name],pdb_resdict[pdb_name],False)
        superimposer.set_atoms(fixed,moving)
        moving_structure=pdb_structures[pdb_name]
        superimposer.apply(moving_structure.get_atoms())
        outfile=open(pdb_name+'_aligned.pdb','w')
        io.set_structure(moving_structure)
        io.save(outfile)
        outfile.close()
        pdb_structures[pdb_name]=moving_structure
    return pdb_structures

def euclidean_distance(atom1_coord,atom2_coord):
    x1,y1,z1=atom1_coord
    x2,y2,z2=atom2_coord
    distance=((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)**0.5
    return distance

def average(numbers):
    return sum(numbers)/float(len(numbers))

def extract_core(best_pdb,pdb_structures,pdb_chains,cutoff):
    io=Bio.PDB.PDBIO()
    best_structure=pdb_structures[best_pdb]
    model=best_structure[0]
    chains=pdb_chains[best_pdb]
    # remove chains not relevant to alignment
    outlying_chains=[]
    for chain in model:
        if not chain.id in chains:
            outlying_chains.append(chain.id)
    for outlying_chain in outlying_chains:
        model.detach_child(outlying_chain)
    # rebuild resdict, necessary because structures have moved due to alignment
    pdb_resdict=make_residues_dict(pdb_structures,pdb_chains)
    best_residues=pdb_resdict[best_pdb]
    # iterate through residues in repr. structure and remove those not, on average, within the cutoff distance
    distances_byresidue=dict()
    for pdb_name in pdb_resdict.keys():
        if pdb_name==best_pdb:
            continue
        # again, have to use alpha carbons... same length not necessary but convenient and
        # might as well recycle function 'get_fixed_moving_CAlists' since it outputs alpha carbons only
        fixed,moving=get_fixed_moving_CAlists(pdb_resdict[best_pdb],pdb_resdict[pdb_name],True)
        for key in fixed.keys():
            atom1_coord=fixed[key].get_coord()
            atom2_coord=moving[key].get_coord()
            # distance calculation, evidently
            distance=euclidean_distance(atom1_coord,atom2_coord)
            try:
                distances=distances_byresidue[key]
            except KeyError:
                distances=[]
            # distances, by residue, across all structure comparison, are stored in a list...
            distances.append(distance)
            # that is then commited to a dictionary where the key is still '*chain*_*position*'
            distances_byresidue[key]=distances
    # loop where residues removed
    for key in distances_byresidue.keys():
        distances=distances_byresidue[key]
        avgdist=average(distances)
        # distance threshold condition
        if avgdist <= cutoff:
            continue
        chainid=key.split('_')[0]
        pos=key.split('_')[1]
        chain=model[chainid]
        res=best_residues[key]
        # residue removal
        chain.detach_child(res.id)
    outfile=open(best_pdb+'_core.pdb','w')
    io.set_structure(best_structure)
    io.save(outfile)
    outfile.close()

def surface_residues_fileparser(surface_residues_filename,chain):
    surface_residues_file=open(surface_residues_filename)
    rows=surface_residues_file.read().strip().split('\n')
    surface_residues=[]
    for row in rows:
        srow=row.strip().split()
        try:
            chan=srow[2].strip()
            pos=srow[3].strip()
        except IndexError:
            continue
        if chain==chan:
            surface_residues.append(pos)
    return surface_residues

def define_contact_residues_fast(structure,chain1id,chain2id,distance_threshold):
    chain1=structure[0][chain1id]
    chain2=structure[0][chain2id]
    chain2ca=[]
    for r in chain2:
        try:
            rca=r['CA']
        except KeyError:
            continue
        chain2ca.append(rca)
    neighbors=Bio.PDB.NeighborSearch(chain2ca)
    contacts1=[]
    contacts2=[]
    for res1 in chain1:
        try:
            r1ca=res1['CA']
        except KeyError:
            continue
        r1ind=res1.get_id()[1]
        r1sym=res1.get_resname()
        for r2ca in neighbors.search(r1ca.get_coord(),distance_threshold):
            res2=r2ca.get_parent()
            r2ind=res2.get_id()[1]
            r2sym=res2.get_resname()
            if not r1ind in contacts1:
                contacts1.append(str(r1ind))
            if not r2ind in contacts2:
                contacts2.append(str(r2ind))
    contacts1=list(set(contacts1))
    contacts2=list(set(contacts2))
    return [contacts1,contacts2]

def define_interface_residues(structure,chainid1,chainid2,distance_threshold):
    contacts1=[]
    contacts2=[]
    chain1=structure[0][chainid1]
    chain2=structure[0][chainid2]
    for res1 in chain1:
        for res2 in chain2:
            found=False
            for atom1 in res1:
                xyz1=atom1.get_coord()
                for atom2 in res2:
                    xyz2=atom2.get_coord()
                    dist=dist_in3d(xyz1,xyz2)
                    if dist<=distance_threshold:
                        contacts1.append(res1.id[1])
                        contacts2.append(res2.id[1])
                        found=True
                        break
                if found:
                    break
    contacts1=list(set(contacts1))
    contacts2=list(set(contacts2))
    return [contacts1,contacts2]

def dist_in3d(xyz1,xyz2):
    x1=xyz1[0]
    y1=xyz1[1]
    z1=xyz1[2]
    x2=xyz2[0]
    y2=xyz2[1]
    z2=xyz2[2]
    dist=((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**(float(1)/2)
    return dist

def define_passive_residues(structure,chainid,contact_residues,surface_residues,distance_threshold,surfaceonly):
    model=structure[0]
    chain=model[chainid]
    contact_CA=[]
    for res in chain:
        respos=str(res.get_id()[1])
        if respos in contact_residues:
            contact_CA.append(res['CA'])
    contact_residues=set(contact_residues)
    try:
        neighbors=Bio.PDB.NeighborSearch(contact_CA)
    except IndexError:
        print contact_CA
        print contact_residues
    passive_residues=[]
    for res in chain:
        try:
            resca=res['CA']
        except KeyError:
            continue
        respos=str(resca.get_parent().get_id()[1])
        for passresca in neighbors.search(resca.get_coord(),distance_threshold):
            passres=passresca.get_parent()
            passrespos=passres.get_id()[1]
            if not respos in contact_residues:
                passive_residues.append(respos)
    if surfaceonly:     
        return list(set(passive_residues).intersection(set(surface_residues)))
    else:
        return list(set(passive_residues))            


#surface_residues_hum=surface_residues_fileparser(docking_dir+'\\'+'1GWP_capsid_passive.rsa','C')
#surface_residues_hiv=surface_residues_fileparser(docking_dir+'\\'+'2CPL_cypa_passive.rsa','G')
#structure=parser.get_structure('1M9C_capsid_cypa',docking_dir+'\\'+'1M9C_capsid_cypa.pdb')
#contacts_humhiv,contacts_hivhum=define_contact_residues_fast(structure,'A','D',6.5)
#passive_residues_hum=define_passive_residues(structure,'A',contacts_humhiv,surface_residues_hum,10.0,True)
#passive_residues_hiv=define_passive_residues(structure,'D',contacts_hivhum,surface_residues_hiv,10.0,True)

def rmsd_ligands(resdict,pdb_structures,pdbname1,pdbname2,chaindict,chainmap,offset):
    resdict1=resdict[pdbname1]
    resdict2=resdict[pdbname2]
    chainmap=dict()
    receptorid1,ligandid1=chaindict[pdbname1]
    receptorid2,ligandid2=chaindict[pdbname2]
    chainmap[receptorid1]=receptorid2
    chainmap[ligandid1]=ligandid2
    fixed,moving=get_fixed_moving_CAlists(resdict1,resdict2,False,chainmap,0)
    #print [res.parent.id[1] for res in fixed]
    #print [res.parent.id[1] for res in moving]
    superimposer=Bio.PDB.Superimposer()
    superimposer.set_atoms(fixed,moving)
    moving_structure=pdb_structures[pdbname2]
    superimposer.apply(moving_structure.get_atoms())
    pdb_structures[pdbname2]=moving_structure
    ligand_chaindict=dict()
    ligand_chaindict[pdbname1]=[ligandid1]
    ligand_chaindict[pdbname2]=[ligandid2]
    resdict=make_residues_dict(pdb_structures,ligand_chaindict)
    resdict1=resdict[pdbname1]
    resdict2=resdict[pdbname2]
    ligand1CA,ligand2CA=get_fixed_moving_CAlists(resdict1,resdict2,False,chainmap,0)
    rmsd_ligands=rmsd_manual(ligand1CA,ligand2CA)
    return rmsd_ligands

def rmsd_interface_residues(resdict,pdb_structures,pdbname1,pdbname2,chaindict,chainmap,offset,ires_r1,ires_l1,ires_r2,ires_l2):
    resdict1=resdict[pdbname1]
    resdict2=resdict[pdbname2]
    chainmap=dict()
    receptorid1,ligandid1=chaindict[pdbname1]
    receptorid2,ligandid2=chaindict[pdbname2]
    chainmap[receptorid1]=receptorid2
    chainmap[ligandid1]=ligandid2
    resdict1_i=resdict_narrow(resdict1,receptorid1,ires_r1,dict())
    resdict1_i=resdict_narrow(resdict1_i,receptorid1,ires_r1,resdict1_i)
    resdict2_i=resdict_narrow(resdict2,receptorid2,ires_r2,dict())
    resdict2_i=resdict_narrow(resdict2_i,receptorid2,ires_r2,resdict2_i)
    fixed,moving=get_fixed_moving_CAlists(resdict1_i,resdict2_i,False,chainmap,0)
    #print fixed
    #print moving
    superimposer=Bio.PDB.Superimposer()
    superimposer.set_atoms(fixed,moving)
    moving_structure=pdb_structures[pdbname2]
    superimposer.apply(moving_structure.get_atoms())
    pdb_structures[pdbname2]=moving_structure
    ichaindict=dict()
    ichaindict[pdbname1]=[resdict1_i]
    ichaindict[pdbname2]=[resdict2_i]
    resdict=make_residues_dict(pdb_structures,chaindict)
    resdict1=resdict[pdbname1]
    resdict2=resdict[pdbname2]
    resdict1_i=resdict_narrow(resdict1,receptorid1,ires_r1,dict())
    resdict1_i=resdict_narrow(resdict1_i,receptorid1,ires_r1,resdict1_i)
    resdict2_i=resdict_narrow(resdict2,receptorid2,ires_r2,dict())
    resdict2_i=resdict_narrow(resdict2_i,receptorid2,ires_r2,resdict2_i)
    ligand1CA,ligand2CA=get_fixed_moving_CAlists(resdict1_i,resdict2_i,False,chainmap,0)
    rmsd_interface=rmsd_manual(ligand1CA,ligand2CA)
    return rmsd_interface

def resdict_narrow(inresdict,chainid,residues,outresdict):
    for key in inresdict.keys():
        chain,pos=key.strip().split('_')
        if chain==chainid and int(pos) in residues:
            outresdict[key]=inresdict[key]
    return outresdict

def percentage_recaptured_interface_residues(ires_r1,ires_l1,ires_r2,ires_l2):
    t1=0
    c1=0
    for ir_r1 in list(ires_r1):
        if ir_r1 in ires_r2:
            c1+=1
        t1+=1
    for ir_l1 in list(ires_l1):
        if ir_l1 in ires_l2:
            c1+=1
        t1+=1
    t2=0
    c2=0
    for ir_r2 in list(ires_r2):
        if not ir_r2 in ires_r1:
            c2+=1
        t2+=1
    for ir_l2 in list(ires_l2):
        if not ir_l2 in ires_l1:
            c2+=1
        t2+=1    
    return [100*float(c1)/t1,100*float(c2)/t2]

##### residuelist1CA and residuelist2CA must be the same length #####
def rmsd_manual(residuelist1CA,residuelist2CA):
    distances=[]
    for i,res1CA in enumerate(residuelist1CA):
        res2CA=residuelist2CA[i]
        x1,y1,z1=res1CA.get_coord()
        x2,y2,z2=res2CA.get_coord()
        distance=((float(x2-x1))**2+(float(y2-y1))**2+(float(z2-z1))**2)**0.5
        distances.append(distance)
    rmsd=(sum([d**2 for d in distances])/(i+1))**0.5                   
    return rmsd

######## Code for scoring HADDOCK runs (HIV recaptured complexes) ########
#pdbname1='1G9N_gp120_CD4'
#pdbname2='CD4_vs_HIV_gp120_cluster1_1'
#chains=[['C','G'],['A','B']]
##pdbname1='3MIA_Tat_PTEFB'
##pdbname2='PTEFB_B_vs_HIV_Tat_cluster2_1'
##chains=[['B','C'],['A','B']]
##pdbname1='1M9C_capsid_cypa'
##pdbname2='CypA_vs_HIV_capsid_cluster2_1'
##chains=[['A','D'],['A','B']]
##pdbname1='1M9C_capsid_cypa'
##pdbname2='CypA_vs_HIV_Capsid_refine_bsa_cluster1_1'
##chains=[['A','D'],['A','B']]
#pdb_filenames=[pdbname1+'.pdb',pdbname2+'.pdb']
#chaindict=dict(zip([pdbname1,pdbname2],chains))
#pdb_structures=import_pdb_structures(pdb_filenames)
#resdict=make_residues_dict(pdb_structures,chaindict)
#chainmap=dict()
#chainmap[chains[0][0]]=chains[1][0]
#chainmap[chains[0][1]]=chains[1][1]
#rmsd_ligands=rmsd_ligands(resdict,pdb_structures,pdbname1,pdbname2,chaindict,chainmap,0)
##### for CD4_gp120 #####
#ires_r1=set([32,33,34,40,41,42,43,44,45,46,47,48,59,60])
#ires_l1=set([124,365,366,367,368,427,429,430,458,459,460,472,473,474])
##### for PTEFB_Tat #####
##ires_r1=set([40,44,47,50,51,53,54,55,56,95,96,97,112,155,172,176,248,249,250])
##ires_l1=set([1,2,3,6,10,11,13,14,15,16,17,36,43,44,46,47,48])
##### for LEDGF_integrase (non-dimer) #####
##ires_r1=set([365,366,407,408])
##ires_l1=set([125,127,128,129,131])
##### for LEDGF_integrase (dimer) #####
##ires_r1=set([364,365,366,367])
##ires_l1=set([1168,1169,1170])
##### for CypA_capsid  #####
##ires_r1=set([102,103,101,73,72])
##ires_l1=set([88,89,86,87,90])
#faceres_r1,faceres_l1=find_interacting_res(pdb_structures[pdbname1][0][chains[0][0]],pdb_structures[pdbname1][0][chains[0][1]],10.0)
#faceres_r2,faceres_l2=find_interacting_res(pdb_structures[pdbname2][0][chains[1][0]],pdb_structures[pdbname2][0][chains[1][1]],10.0)
#faceres_r1=set(faceres_r1)
#faceres_l1=set(faceres_l1)
#faceres_r2=set(faceres_r2)
#faceres_l2=set(faceres_l2)
#ires_r2,ires_l2=define_contact_residues_fast(pdb_structures[pdbname2],'A','B',6.5)
#ires_r2=[int(pos) for pos in ires_r2]
#ires_l2=[int(pos) for pos in ires_l2]
#ires_r2=set(ires_r2)
#ires_l2=set(ires_l2)
#rmsd_interface=rmsd_interface_residues(resdict,pdb_structures,pdbname1,pdbname2,chaindict,chainmap,0,faceres_r1,faceres_l1,faceres_r2,faceres_l2)
#recapture=percentage_recaptured_interface_residues(ires_r1,ires_l1,ires_r2,ires_l2)

##### Code to calculate interface residues for a protein complex #####
#os.chdir('E:\\SNPs@interface\\Docking\\Mutagenesis_docking\\structures\\beta_trypsin_BPTI')
#pdb_filenames=['beta_trypsin_BPTI.pdb']
#pdb_structures=import_pdb_structures(pdb_filenames)
#contacts_enzyme,contacts_ligand=find_interacting_res(pdb_structures['beta_trypsin_BPTI'][0]['E'],pdb_structures['beta_trypsin_BPTI'][0]['I'],10.0)
