import re
import numpy as np
import pandas as pd
from ..coords.base import CoordsIO
from .script import Script
from ..utils.strings import clean
from ..utils.errors import ParserError, MiMiCPyError
from ..utils.constants import BOHR_RADIUS

class Pseudopotential:

    def __init__(self, coords, pp_type='MT_BLYP.psp', labels='KLEINMAN-BYLANDER', lmax='S', loc=''):
        if all(isinstance(i, list) for i in coords):
            self.coords = coords
        else:
            self.coords = [coords]
            
        self.pp_type = pp_type
        self.labels = labels
        self.lmax = lmax
        self.loc = loc

    def __str__(self):
        if not self.pp_type.startswith('_'):
            self.pp_type = '_' + self.pp_type
        if not self.labels.startswith(' ') and self.labels != '':
            self.labels = ' ' + self.labels

        pp_block = '{}{}\n'.format(self.pp_type, self.labels)
        pp_block += '    LMAX={}'.format(self.lmax.upper())
        pp_block += '\n' if self.loc == '' else ' LOC={}\n'.format(self.loc.upper())
        pp_block += '    {}\n'.format(len(self.coords))

        for row in self.coords:
            pp_block += ' {:>18.12f} {:>18.12f} {:>18.12f}\n'.format(row[0], row[1], row[2])

        return pp_block

    @classmethod
    def from_string(cls, text):
        splt = text.splitlines()
        
        pp, labels = splt[0].split(' ', 1)
        labels = labels.strip()
        
        second_line_regex = re.compile(r'(\w+)\s*=\s*(\w)')
        second_line = dict(second_line_regex.findall(splt[1].upper()))
        
        if 'LMAX' in second_line:
            lmax = second_line['LMAX']
        else:
            raise ParserError(file='CPMD script', details='no Lmax data for atom')
        
        if 'LOC' in second_line:
            loc = second_line['LOC']
        else:
            loc = ''
        
        no = int(splt[2].strip())
        string_to_list = lambda string: [float(i) for i in string.split()]
        coords = [string_to_list(i) for i in splt[3:]]
        
        if len(coords) != no:
            raise ParserError(file='CPMD script', details='mismatch in no. of atoms ({} vs {})'.format(no, len(coords)))     
            
        return cls(coords, pp, labels, lmax, loc)


class Section(Script):

    def __str__(self):
        section_string = ''
        for keyword in self.parameters:
            value = getattr(self, keyword)
            if isinstance(value, Pseudopotential):
                section_string += '\n*{}{}'.format(keyword.title(), str(value))
            else:
                section_string += '\n    {}'.format(keyword.upper().replace('__', ' ').replace('_', '-'))
                if keyword == 'overlaps':
                    value = '\n'.join(['        '+s.strip() for s in value.splitlines()])
                    section_string += '\n{}'.format(str(value))
                elif value is not True:
                    section_string += '\n        {}'.format(str(value))
            
        return section_string
    
    @staticmethod
    def __chknumeric(s):
        splt = s.split()
        
        if len(splt) == 1:
            return s.replace('.','').replace('-','').isnumeric()
        else:
            for i in splt:
                if not Section.__chknumeric(i):
                    return False
            return True

    @classmethod
    def from_string(cls, text, as_atoms=False):
        
        if as_atoms: return Section.__from_string_atoms(text)
        
        i = 0
        section = cls()
        splt = text.splitlines()
        
        if len(splt) == 1: setattr(section, splt[i], True)
        
        while i < len(splt):
            splt_i = splt[i].strip()
            
            if splt_i == 'PATHS':
                try:
                    setattr(section, splt_i, "\n".join([s.strip() for s in splt[i+1:i+3]]))
                except IndexError:
                    raise ParserError(file='CPMD script', details='PATHS in MIMIC section not formatted correctly')
                i += 2
                
            elif splt_i == 'OVERLAPS':
                try:
                    no = int(splt[i+1])
                    ov = "\n".join([s.strip() for s in splt[i+1:i+no+2]])
                except IndexError:
                    raise ParserError(file='CPMD script', details='OVERLAPS in MIMIC section not formatted correctly')
                setattr(section, splt_i, ov)
                
                i += no+1  
                    
            elif i < len(splt)-1 and Section.__chknumeric(splt[i+1].strip()):            
                setattr(section, splt_i, splt[i+1].strip())
                i += 1
                
            elif not Section.__chknumeric(splt[i]):
                setattr(section, splt_i, True)
                
            i += 1
            
        return section
    
    @classmethod
    def __from_string_atoms(cls, text):
        section = cls()
        
        ids = [ (i.start(), i.end()) for i in re.finditer('\*', text)]
        
        for i, _ in enumerate(ids):
            if i == len(ids)-1:
                next_id = -1
            else:
                next_id = ids[i+1][0]
                
            atom_txt = text[ids[i][1]:next_id]
            elem, a_txt = atom_txt.split('_', 1)
            setattr(section, elem, Pseudopotential.from_string(a_txt))
        
        return section


class CpmdScript(Script):

    def __init__(self, *sections):
        super().__init__()
        for section in sections:
            setattr(self, section, Section())

    def __str__(self):
        cpmd_script = ''
        for section in self.parameters:
            if section == 'info':
                section_string = '\n'+getattr(self, section)
            else:
                section_string = str(getattr(self, section))
            cpmd_script += '\n&{}{}\n&END\n'.format(section.upper(), section_string)
        return cpmd_script
    
    @classmethod
    def from_string(cls, text):
         text = clean(text)
         section_reg = re.compile(r'\s*\&(.*?)\n((?:.+\n)+?)\s*(?:\&END)')
         sections = section_reg.findall(text)
         
         inp = cls()
         
         for k,v in sections:
             if k == 'INFO': setattr(inp, k, v.strip())
             elif k == 'ATOMS': setattr(inp, k, Section.from_string(v, as_atoms=True))
             else: setattr(inp, k, Section.from_string(v))
         
         
         return inp
     
    def to_coords(self, mpt, out, title=None, ext=None):
        if not self.has_parameter('mimic'):
            raise MiMiCPyError('MIMIC section not found in CPMD script')
        elif not self.mimic.has_parameter('overlaps'):
            raise MiMiCPyError('OVERLAPS in MIMIC section not found in CPMD script')
            
        try:
            ids = [int(i.split()[1]) for i in self.mimic.overlaps.splitlines()[1:]]
        except (ValueError, IndexError):
            raise ParserError(file='CPMD script', details='OVERLAPS in MIMIC section not formatted correctly')
        
        if not self.has_parameter('atoms'):
            raise MiMiCPyError('No ATOMS section found in CPMD script')
        
        coords_list = []
        for k,v in self.atoms.parameters.items():
            coords_list  += v.coords
        
        coords_np = np.array(coords_list)*BOHR_RADIUS
        
        if len(ids) != coords_np.shape[0]:
            raise MiMiCPyError('Mismatch between no. of atoms in OVERLAPS and ATOMS sections ({} vs {})'.format(len(ids), coords_np.shape[0]))
        
        coords = pd.DataFrame({'id': ids, 'x': coords_np[:,0], 'y': coords_np[:,1], 'z': coords_np[:,2]})
        
        if not title: title = 'Coordinates from CPMD/MiMiC script'
        
        CoordsIO(out, mode='w', ext=ext).write(mpt, coords, title=title)