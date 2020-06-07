import pandas as pd
import re
from ..utils.constants import element_names
from .._global import _Global as gbl
from .parser import Parser
from ..utils.errors import ParserError
 
def isSection(section, txt):
    if section == '*': section = ''
    if '[' in txt and ']' in txt and section in txt:
        return True
    else: return False
    
def getSection(section, txt):
    # txt is assumed to be clean
    # i.e., no comments or double new lines
    
    # find text b/w [ section ] and either [ or # (for #if, etc.) or EOF
    reg = re.compile(fr"\[\s*{section}\s*\]\n((?:.+\n)+?)(?:$|\[|#)", re.MULTILINE)
    r = reg.findall(txt)
    return r

def cleanText(txt):
    txt_ = re.sub(re.compile(";(.*)\n" ) ,"\n" , txt) # strip comments
    return re.sub(re.compile("\n{2,}" ) ,"\n" , txt_) # remove double new lines
        
def molecules(tail):    
    _mols = []
    
    for line in tail.splitlines()[::-1]: # traverse backwards
        if isSection('molecules', line):
            break
        elif line.strip() == '' or line.startswith(';'):
            continue
        mol, no = line.split()
        _mols += [(mol, int(no))]

    return _mols[::-1]

def _read_atomtypes(itp_file, buff):
    """ Function to read file for atomtypes by chunks"""
    atomtypes = ''
    end = ['bondtypes', 'moleculetypes']
    
    file = Parser(itp_file, buff)
    for chunk in file:
        atomtypes += chunk
        if any([isSection(hdr, chunk) for hdr in end]): break
    
    atomtypes = cleanText(atomtypes)
    return getSection('atomtypes', atomtypes)[0]

def atomtypes(itp_file, buff):
    
    atomtypes = _read_atomtypes(itp_file, buff)
    
    atm_types_to_symb = {} # init atom types to symbol
    
    for line in atomtypes.splitlines():
        e = line.split()[0] # first val is atom type
        _n = line.split()[1]
        
        if not _n.isnumeric():
            # should raise an exception here
            continue
        else: n = int(_n)-1 # second is element no.
        
        if n == -1: continue # dummy masses, skip for now
        atm_types_to_symb[e]  = element_names[n] # fill up at
    
    return atm_types_to_symb

def non_std_atomtypes(itp_file, buff):
    atomtypes = _read_atomtypes(itp_file, buff)
    return [line.split()[0] for line in atomtypes.splitlines()]


class ITPParser:    
    
    columns = ['number', 'type', 'resid','resname','name', 'charge','element',	'mass']
    dfs = []
    mols = []

    def __init__(self, mols_to_read, atm_types_to_symb, buff, guess):
        self.mols_to_read = mols_to_read
        self.atm_types_to_symb = atm_types_to_symb
        self.buff = buff
        self.guess = guess
        
    def read(self, itp_file):
        
        file = Parser(itp_file, self.buff)
        
        molecule = False
        txt = ''
        
        for chunk in file:
            if isSection('moleculetype', chunk):
                molecule = True
                
            if molecule: txt += chunk
            
            if isSection('bonds', chunk):
                molecule = False
        
        return txt
    
    def parse(self, file_name, itp_text=None):
        
        if itp_text == None:
            itp_text = self.read(file_name)
        
        itp_text = cleanText(itp_text)
        
        mol_section = getSection('moleculetype', itp_text)
        atom_section = getSection('atoms', itp_text)
        
        for m, a in zip(mol_section, atom_section):
            mol = m.split()[0]
            if mol not in self.mols_to_read: continue
            self.mols.append(mol)
            self._parseatoms(a, file_name)
        
    def _parseatoms(self, txt, file_name):
        
        df_ = {k:[] for k in self.columns}
        
        for line in txt.splitlines():
            
            splt = line.split()
            if len(splt) == 8:
                nr, _type, resnr, res, name, cgnr, q, mass = splt[:8]
            elif len(splt) == 7:
                nr, _type, resnr, res, name, cgnr, q = splt[:7]
                mass = 0
            else:
                continue
            
            c = self.columns
            df_[c[0]].append(int(nr))
            df_[c[1]].append(_type)
            df_[c[2]].append(int(resnr))
            df_[c[3]].append(res)
            df_[c[4]].append(name)
            df_[c[5]].append(float(q))
            
            if _type in self.atm_types_to_symb:
                elem = self.atm_types_to_symb[_type]
            elif self.guess:
                mass_int = int(float(mass))
                
                if mass_int <= 0:
                    raise ParserError(file=file_name, ftype="topolgy", \
                                     extra=(f"Cannot determine atomic symbol for atom ID {nr} and name {name} as mass"
                                             "information is not available from the force field"))
                # guess atomic no from mass
                # works well if no isotopes present
                
                if mass_int<=1: elem = 'H' # for H
                elif mass_int<36: elem = element_names[mass_int//2 - 1] # He to Cl
                else: elem = name.title() # from Ar onwards, assume name same as symbol, case insensitive
                
                gbl.logger.write('warning', (f"Guessing atomic symbol for atom id {nr} and name {name} as {elem}..") )
            
            else:
                raise ParserError(file=file_name, ftype="topology", \
                                     extra="Cannot determine atomic symbol for atom ID {nr} and name {name}")
            
            df_[c[6]].append(elem)
            df_[c[7]].append(mass)
        
        df = pd.DataFrame(df_).set_index(['number'])
        self.dfs.append( [ len(df), df ] )