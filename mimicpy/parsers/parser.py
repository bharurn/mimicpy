from .._global import _Global as gbl

class Parser(object):
    def __init__(self, file, buff):
        self.f = gbl.host.open(file, 'rb')
        self.buff = buff
        
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.next()

    def next(self):
        out = self.f.read(self.buff)
    
        if out == b'':
            self.close()
            raise StopIteration()
        else: return out.decode()
    
    def __del__(self):
        self.close()
    
    def close(self):
        self.f.close()