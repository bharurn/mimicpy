from collections import OrderedDict 

class Section:
    def __init__(self, **kwargs):
        
        for key, value in kwargs.items():
            setattr(self, key, value)
        
    def __str__(self):
        data = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith('__')]
        
        val = ''
        
        for d in data:
            if getattr(self, d) == None:
                continue
            
            d_ = d.replace('_','-').replace('--', ' ').upper()
            v = str(getattr(self, d)).strip()
            if v == '':
                val += d_+'\n'
            else:
                val += f"{d_}\n{v}\n"
        
        return val
    
class Atom:
    def __init__(self, coords=[], lmax='s', pp='MT_BLYP', labels=''):
        self.coords = coords
        self.pp = pp
        self.labels = labels
        self.lmax = lmax
    
    def __str__(self):
        if not self.pp.startswith('_'): self.pp = '_' + self.pp
        if not self.labels.startswith(' ') and self.labels != '': self.labels = ' ' + self.labels
        val = f'{self.pp}{self.labels}\n'
        val += f'   LMAX={self.lmax.upper()}\n'
        val += f'    {len(self.coords)}\n'
        for d in self.coords:
            val += f'  {d[0]}   {d[1]}   {d[2]}\n'
        val += '\n'
        return val
    
class Input:
    def __init__(self, *args):
        self.atoms = OrderedDict()
        self.info = 'MiMiC Run'
        for val in args:
            setattr(self, val, Section())
    
    def checkSection(self, section):
        if callable(section) or section.startswith('__'): return False
        
        return hasattr(self, section)
    
    def __str__(self):
        data = [attr for attr in dir(self) if not callable(getattr(self, attr)) and not attr.startswith('__')]
        
        val = ''
        
        for d in data:
            if getattr(self, d) == None:
                continue
            elif d.upper() == 'INFO':
                val += f'\n&INFO\n{getattr(self, d)}\n&END'
            elif d.upper() == 'ATOMS':
                val += '\n&ATOMS\n'
                for k, v in getattr(self, d).items():
                    val += f'*{k}{str(v)}'
                val += '&END\n'
            else:
                v = str(getattr(self, d))
                val += f"\n&{d.upper()}\n{v}&END\n"
        
        return val