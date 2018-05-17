import os
import sys
os.environ['PYMOL_PATH']='C:\\Program Files (x86)\\PyMOL\\PyMOL\\PyMOL.exe'
# import pymol command line interface module, cmd
import cmd
# change working directory to one containing structral files to be analyzed
os.chdir('C:\\Users\\Norman\\Documents\\PARABON')
# add the directory for Python packages to system path, so that import function can find Bio.PDB
sys.path.append('C:\\Python27\\Lib\\site-packages')
# import Bio.PDB.PDBParser, allows non-pymol analyses to be performed
import Bio.PDB.PDBParser

class protein_complex:
    def __init__(self,complexname,receptor,ligand):
        try:
            cmd.load(complexname)
        except CmdException:
            print complexname
        self.name=complexname.split('.')[0]
        self.chains=dict()
        try:
            chainids=cmd.get_chains(self.name)
        except CmdException:
            print self.name
        for chainid in chainids:
            self.chains[chainid]=cmd.select(self.name+'_'+chainid,self.name+' and chain '+chainid)
            #cmd.delete(self.name+'_'+chainid)
        self.receptor=receptor
        self.ligand=ligand
    def find_interfaces(self):
        interfaces=dict()
        chainkeys=self.chains.keys()
        for chain1key in chainkeys:
            chain1=self.chains[chain1key]
            for chain2key in chainkeys:
                if chain1key==chain2key:
                    continue
                chain2=self.chains[chain2key]
                interface1=cmd.select(self.name+'_'+chain1key+'_interface_'+chain1key+chain2key,
                                      'byres '+self.name+' and chain '+chain1key+' within 4.0 of chain '+chain2key)
                interface2=cmd.select(self.name+'_'+chain2key+'_interface_'+chain1key+chain2key,
                                      'byres '+self.name+' and chain '+chain2key+' within 4.0 of chain '+chain1key)
                newinterface=cmd.select(self.name+'_interface_'+chain1key+chain2key,
                                                          '('+self.name+'_'+chain1key+'_interface_'+chain1key+chain2key+','
                                                          +self.name+'_'+chain2key+'_interface_'+chain1key+chain2key+')')
                interfaces[chain1key+chain2key]=newinterface
                cmd.delete(self.name+'_'+chain1key+'_interface_'+chain1key+chain2key)
                cmd.delete(self.name+'_'+chain2key+'_interface_'+chain1key+chain2key)
        self.interfaces=interfaces

class compare_complexes(protein_complex):
    def __init__(self,complex_names):
        protein_complexes=dict()
        for name in complex_names.keys():
            receptor,ligand=complex_names[name]
            newcomplex=protein_complex(name,receptor,ligand)
            shortname=name.split('.')[0]
            protein_complexes[shortname]=newcomplex
        self.complexes=protein_complexes
    def rmsd_interfaces(self):
        mat=[]
        for i,complex1name in enumerate(self.complexes.keys()):
            complex1=self.complexes[complex1name]
            r1=complex1.receptor
            l1=complex1.ligand
            row=[]
            for j,complex2name in enumerate(self.complexes.keys()):
                if i==j:
                    row.append(float(0))
                    continue
                complex2=self.complexes[complex2name]
                r2=complex2.receptor
                l2=complex2.ligand
                complex1.find_interfaces()
                complex2.find_interfaces()
                interface1=complex1.name+'_interface_'+r1+l1
                interface2=complex2.name+'_interface_'+r2+l2
                rmsd=float(cmd.align(interface1,interface2)[0])
                row.append(rmsd)
            mat.append(row)
        self.rmsd_interfaces=mat
    def rmsd_ligands(self):
        mat=[]
        for i,complex1name in enumerate(self.complexes.keys()):
            complex1=self.complexes[complex1name]
            r1=complex1.receptor
            l1=complex1.ligand
            row=[]
            for j,complex2name in enumerate(self.complexes.keys()):
                if i==j:
                    row.append(float(0))
                    continue
                complex2=self.complexes[complex2name]
                r2=complex2.receptor
                l2=complex2.ligand
                rmsd=cmd.align(complex1name+' and chain '+l1,complex2name+' and chain '+l2)[0]
                row.append(rmsd)
            mat.append(row)
        self.rmsd_ligands=mat
    def rmsd_ligands_ralign(self):
        parser=Bio.PDB.PDBParser()
        mat=[]
        for i,complex1name in enumerate(self.complexes.keys()):
            print complex1name
            complex1=self.complexes[complex1name]
            r1=complex1.receptor
            l1=complex1.ligand
            row=[]
            for j,complex2name in enumerate(self.complexes.keys()):
                print complex2name
                if i==j:
                    row.append(float(0))
                    continue
                complex2=self.complexes[complex2name]
                r2=complex2.receptor
                l2=complex2.ligand
                cmd.align(complex1name+' and chain '+r1,complex2name+' and chain '+r2)
                outobj1=complex1name.split('_')[0]+'_ligand'
                cmd.create(outobj1,complex1name+' and chain '+l1)
                cmd.save(outobj1+'.pdb',outobj1)
                outobj2=complex2name.split('_')[0]+'_ligand'
                cmd.create(outobj2,complex2name+' and chain '+l2)
                cmd.save(outobj2+'.pdb',outobj2)
                structure1=parser.get_structure(outobj1,outobj1+'.pdb')
                model1=structure1[0]
                ligand1=model1[l1]
                structure2=parser.get_structure(outobj2,outobj2+'.pdb')
                model2=structure2[0]
                ligand2=model2[l2]
                distances=[]
                l1_residues=dict()
                l2_residues=dict()
                for res1 in ligand1:
                    pos1=str(res1.get_id()[1])
                    l1_residues[pos1]=res1
                for res2 in ligand2:
                    pos2=str(res2.get_id()[1])
                    l2_residues[pos2]=res2
                for pos1 in l1_residues.keys():
                    try:
                        res2=l2_residues[pos1]
                    except KeyError:
                        continue
                    res1=l1_residues[pos1]
                    try:
                        x1,y1,z1=res1['CA'].get_coord()
                        x2,y2,z2=res2['CA'].get_coord()
                    except KeyError:
                        continue
                    distance=((float(x2-x1))**2+(float(y2-y1))**2+(float(z2-z1))**2)**0.5
                    distances.append(distance)
                rmsd=sum([d**2 for d in distances])**0.5
                row.append(rmsd)
            mat.append(row)
        self.rmsd_ligands_ralign=mat
        print "done"
        #delete file with aligned receptors

pdb_filenames=['C:\\Users\\Norman.Goodacre\\Documents\\Ecoli_Hpr_docking\\HADDOCK\\3CCD.fixed.A_4YNG.fixed.A_cluster1_1.pdb','C:\\Users\\Norman.Goodacre\\Documents\\Ecoli_Hpr_docking\\ZDOCK\\job_111614_3CCD.fixed.A_4YNG.fixed.A\\complex.1.pdb']
chains=[['A','B'],['B','A']]
complexes=compare_complexes(dict(zip(pdb_filenames,chains)))
complexes.rmsd_interfaces()
complexes.rmsd_ligands()
complexes.rmsd_ligands_ralign()


