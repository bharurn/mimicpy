from ..system import _hndlpdb as hpdb
from ..utils.scripts import mdp
from .base import BaseHandle
from .simulate import MD
import mimicpy._global as _global
from . import _qmhelper, _mpt_helper
from ._constants import hartree_to_ps
from ..scripts import cpmd
from collections import defaultdict
import pandas as pd
import pickle

class MM(BaseHandle):
    
    def __init__(self, status=[]):
        self._ion_kwargs = {'pname': 'NA', 'nname': 'CL', 'neutral': ''}
        self._box_kwargs = {'c': '', 'd': 1.9, 'bt': 'cubic'}
        self._solavte_kwargs = {'cs': 'spc216.gro'}
        self._topol_kwargs = {'water': 'tip3p', 'ff': 'amber99sb-ildn'}
        self.his_str = ''
        super().__init__(status)
        self.dir = 'prepareMM'
        
        self.confin = "confin.pdb"
        self.conf = 'conf.pdb'
        self.conf1 = 'conf1.gro'
        self.conf2 = 'conf2.gro'
        self.conf3 = 'conf3.gro'
        self.topol = f"topol.top"
        self.ions_mdp = f"ions.mdp"
        self.ions_tpr = f"ions.tpr"
    
    def topolParams(self, **kwargs):
        if 'his' in kwargs:
            self.his_str = '\n'.join(kwargs['his'])
            self.his_str = self.his_str.replace('D1', '0')
            self.his_str = self.his_str.replace('E2', '1')
            kwargs['his'] = ''
        
        self._topol_kwargs = kwargs
        
    def getTopology(self, protein):
        self.setcurrent(key='prepMM')
        
        print('Preparing protein topology..')
        
        print(f"Writing {protein.name} to confin.pdb..")
        
        _global.host.write(protein.pdb+protein.water, f'{self.dir}/{self.confin}')
        
        if self.his_str == '':
            print("Letting Gromacs calculate histidine protonantion states..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs)
        else:
            print("Reading histidine protonation states from list..")
            self.gmx('pdb2gmx', f = self.confin, o = self.conf, dirc=self.dir, **self._topol_kwargs, stdin=self.his_str)
        
        conf = f'{self.dir}/{self.conf}'
        
        pdb = _global.host.read(conf)
        
        lines = []
        df_lst = []
        splt = pdb.splitlines()
        for i, line in enumerate(splt[::-1]):
            vals = hpdb.readLine(line)
            if vals['record'] != 'HETATM' and vals['record'] != 'ATOM':
                lines.append(line) 
            elif vals['record'] == 'HETATM' or vals['record'] == 'ATOM':
                if vals['resName'] == 'HOH': lines.append(line)
                else:
                    #create topology dataframe
                    vals.update({'number':i+1})
                    df_lst.append(vals)
        
        df = pd.DataFrame(df_lst)
        df = df.drop(['serial','record','altLoc','resSeq','iCode','occupancy','x','y','z','tempFactor','charge'], axis=1)
        df = df.set_index(['number'])
        
        print("Saving topology data as pickled dataframe..")
        f = _global.host.vi(f'{self.dir}/topol.mpt', 'wb')
        pickle.dump(df, f)

        print("Adding non-standard residues to structure..")
        conf_pdb = '\n'.join(splt[:len(splt)-i]) + '\n' + protein.ligand_pdb + '\n'.join(lines[::-1])
        
        _global.host.write(conf_pdb, conf)
        
        top1 = (f"; Topology data for all non-standard resiudes in {protein.name} created by MiMiCPy\n"
            "; AmberTools was used to generate topolgy parameter for Amber Force Field, conversion to GMX done using Acpype"
                    "\n\n[ atomtypes ]\n")
        top2 = ''
        
        print("Combining ligands topology into single file ligands.itp..")
        for ligname, lig in protein.ligands.items():
            t1, t2 = lig.splitItp()
	
            top1 += t1
            top2 += t2 
            
            print(f"Writing position restraint file for {lig.name}")
            
            _global.host.write(lig.posre, f"{self.dir}/posre_{lig.name}.itp")
            
            _global.host.run(f'echo {lig.name} {lig.chains} >> topol.top')
        
        _global.host.write(top1+'\n'+top2, f"{self.dir}/ligands.itp")
        
        topol = f"{self.dir}/{self.topol}"
        
        _global.host.run(r'sed -i -r "/^#include \".+.ff\/forcefield.itp\"/a #include \"ligands.itp\"" '+topol)
        _global.host.run('grep -v SOL topol.top > topol_.top && mv topol_.top '+topol)
        _global.host.run(f'echo SOL {protein.hoh_mols} >> '+topol)
        print("ligands.itp added to topol.top")
        
        print('Topology prepared..')
        
        self.saveToYaml()
        
    def ionParams(self, **kwargs): self._ion_kwargs = kwargs
    def boxParams(self, **kwargs): self._box_kwargs = kwargs
    def solvateParams(self, **kwargs): self._solavte_kwargs = kwargs
    
    def getBox(self, genion_mdp = mdp.MDP.defaultGenion()):
        print("Preparing system box..")
        print("Solavting protein..")
        
        self.gmx('editconf', f = self.conf, o = self.conf1, dirc=self.dir, **self._box_kwargs)
        
        self.gmx('solvate', cp = self.conf1, o = self.conf2, p = self.topol, dirc=self.dir, **self._solavte_kwargs)
        
        _global.host.write(str(genion_mdp), f'{self.dir}/{self.ions_mdp}')
        
        print("Adding ions to neutralize charge..")
        self.gmx('grompp', f = self.ions_mdp, c = self.conf2, p = self.topol, o = self.ions_tpr, dirc=self.dir)
        
        self.gmx('genion', s = self.ions_tpr, o = self.conf3, p = self.topol, dirc=self.dir, **self._ion_kwargs, stdin="SOL")
        
        print('Simulation box prepared..')
        
        self.saveToYaml()
        
class QM(MD):
    
    def __init__(self, status=defaultdict(list)):
        super().__init__(status)
        
        mpt = self.getcurrent('mpt')
        print(f"Reading topology from {mpt}..")
        df_f = _global.host.vi(mpt, 'rb')
        self.df = pickle.load(df_f)
        
        coords = self.getcurrent('gro') # TO DO: check if latest run is trr or gro, and if trr convert
        print(f"Combining with latest coordinates data from {coords}..")
        self.df, self._mm_box = _mpt_helper.combine(self.df, coords)
        
        self.inp = cpmd.Input()
        
        self.qmatoms = None
        
        self.dir = 'prepareQM'
        
        self.index = 'index.ndx'
        self.preprc = 'processed.top'
        self.mdp_tpr = 'mimic'
        
    def add(self, selection, link=False):
        if isinstance(selection, str):
            qdf = _qmhelper.parse_selec(selection, self.df)
        else:
            qdf = self.df[selection(self.df)]
        qdf.insert(2, 'link', [int(link)]*len(qdf))
        if self.qmatoms is None:
            self.qmatoms = qdf
        else:
            self.qmatoms = self.qmatoms.append(qdf)
    
    def getInp(self, mdp, inp=cpmd.Input()):
        dirc = self.dir
        self.setcurrent(key='prepQM')
        
        print("Changing Gromacs integrator to MiMiC..")
        #mdp.integrator = 'mimic'
        
        print(f"Writing atoms in QM region to {self.index}..")
        #mdp.QMMM_grps = 'QMatoms'
        _global.host.write(_qmhelper.index(self.qmatoms.index), f'{dirc}/{self.index}')
        
        print("Generating Gromacs tpr file for MiMiC run..")
        self.grompp(mdp, f'{self.mdp_tpr}', pp=self.preprc, n=self.index, dirc=dirc)
        
        unk_lst = self.qmatoms[self.qmatoms['element'].apply(lambda x: False if x.strip() != '' else True)]['name'].to_list()
        
        print("Reading force field data to fill in charges and missing atomic symbol information..")
        # the full processed.top file is loaded into memory for fast manipulation, may not be feasible for large files
        conv, q_conv = _qmhelper.pptop(unk_lst, self.qmatoms['name'].to_list(), \
                                       _global.host.read(f'{dirc}/{self.preprc}') ) # read element conv and charge conv
       
        # iterate through df to combine name & element, get charges for each name also
        elems, q = zip(*self.qmatoms[['name', 'element']].apply(\
                       lambda x: (x[1], q_conv[x[0]]) if x[1].strip() != '' else (conv[x[0]], q_conv[x[0]]), axis=1))

        qm_df = self.qmatoms.assign(element=elems)
        self.qmatoms = qm_df # reassign it to original df
        
        self.qmatoms = self.qmatoms.sort_values(by=['link', 'element']).reset_index()
        
        print("Creating CPMD input script..")
        inp.mimic = cpmd.Section()
        inp.mimic.paths = "1\n---" #path will be set in MiMiC run function
        inp.mimic.box = '  '.join([str(s) for s in self._mm_box])
        
        inp = _qmhelper.getOverlaps_Atoms(self.qmatoms, inp)
        
        inp.system.charge = round(sum(q), 1) # system section already created in getQverlap_Atoms()
        
        if not inp.checkSection('cpmd'): inp.cpmd = cpmd.Section()
        inp.cpmd.mimic = ''
        inp.cpmd.parallel__constraints = ''
        
        # set no of steps from mdp file
        if not mdp.hasparam('nsteps'): mdp.nsteps = 1000
        inp.cpmd.maxsteps = mdp.nsteps
        
        # set timestep from mdp file
        if not mdp.hasparam('dt'): mdp.dt = 0.0001
        inp.cpmd.timestep = round(mdp.dt/hartree_to_ps)
        
        self.inp = inp
        
        self.toYaml()
        
        return inp