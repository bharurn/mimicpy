import pandas as pd
import re
from ..utils.constants import element_names
from .._global import _Global as gbl
from .parser import Parser
from ..utils.errors import ParserError

include_file_regex = re.compile(r"#include\s+[\"\'](.+)\s*[\"\']", re.MULTILINE)
 
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

def getSectionHeaders(txt):
    reg = re.compile(fr"\[\s*(.*?)\s*\]", re.MULTILINE)
    r = reg.findall(txt)
    return r

def parseBlocktillSection(file, *sections):
    #parse file till sections (its a list of sections) is found
    itp_txt = ''
    sections_read = [] #keep track of sections read
    
    for chunk in file:
        itp_txt += chunk
        
        sections_read += getSectionHeaders(chunk)
        
        # if all sections have been read, and we now have read another sections after this
        # means that we read all sections we wanted, so break
        if set(sections).issubset(sections_read) and sections_read[-1] != sections[-1]:
            break
        
    return itp_txt

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

def atomtypes(itp_file, buff):
    
    file = Parser(itp_file, buff)
    itp_txt = cleanText(parseBlocktillSection(file, 'atomtypes'))
    
    atomtypes_txt = getSection('atomtypes', itp_txt)
    
    if atomtypes_txt == []:
        # if no atomtypes section look for #include
        for i in include_file_regex.findall(itp_txt):
            full_path = gbl.host.join(gbl.host.dirname(itp_file), i)
            atm_types_to_symb = atomtypes(full_path, buff) # recursively call this func
            if atm_types_to_symb != {}: return atm_types_to_symb # if found atomtypes, return it
        return {} # if not found return empty dict
    else:
        atomtypes_txt = atomtypes_txt[0]
    
    atm_types_to_symb = {} # init atom types to symbol
    
    for line in atomtypes_txt.splitlines():
        e = line.split()[0] # first val is atom type
        _n = line.split()[1]
         
        if not _n.isnumeric():
            # should raise an exception here
            continue
        else: n = int(_n)-1 # second is element no.
        
        if n == -1: continue # dummy masses, skip for now
        atm_types_to_symb[e]  = element_names[n] # fill up at
    
    return atm_types_to_symb

class ITPParser:    
    
    columns = ['number', 'type', 'resid', 'resname', 'name', 'charge', 'element', 'mass']
    dfs = []
    mols = []

    @staticmethod
    def clear():
        ITPParser.mols = []
        ITPParser.dfs = []

    def __init__(self, mols_to_read, atm_types_to_symb, buff, guess):
        self.mols_to_read = mols_to_read
        self.atm_types_to_symb = atm_types_to_symb
        self.buff = buff
        self.guess = guess
    
    def read(self, file_name):
        file = Parser(file_name, self.buff)
        
        itp_text = ''
        
        while not file.isclosed:
            # parse all moleculetype/atoms sections
            chunk = parseBlocktillSection(file, 'moleculetype', 'atoms')
            
            if any([isSection(i, chunk) for i in ['moleculetype', 'atoms']]):
                itp_text += chunk
            
        return itp_text
    
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
            mass = float(mass)
            
            if _type in self.atm_types_to_symb:
                elem = self.atm_types_to_symb[_type]
            elif self.guess:
                mass_int = int(mass)
                
                if mass_int <= 0:
                    raise ParserError(file=file_name, ftype="topolgy", \
                                     extra=(f"Cannot determine atomic symbol for atom ID {nr} and name {name} in residue {res} as mass"
                                             "information is not available from the force field"))
                # guess atomic no from mass
                # works well if no isotopes present
                
                if mass_int<=1: elem = 'H' # for H
                elif mass_int<36: elem = element_names[mass_int//2 - 1] # He to Cl
                else: elem = name.title() # from Ar onwards, assume name same as symbol, case insensitive
                
                gbl.logger.write('warning', (f"Guessing atomic symbol for atom id {nr} and name {name} in residue {res} as {elem}..") )
            else:
                raise ParserError(file=file_name, ftype="topology", \
                                     extra=f"Cannot determine atomic symbol for atom ID {nr} and name {name} in residue {res}")
            
            df_[c[6]].append(elem)
            df_[c[7]].append(mass)
        
        df = pd.DataFrame(df_).set_index(c[0])
        self.dfs.append(df)