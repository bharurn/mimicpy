from .._global import _Global as gbl

class Parser(object):
    def __init__(self, file, buff, decode=True):
        self.f = gbl.host.open(file, 'rb')
        self.buff = buff
        self.decode = decode
        self.isclosed = False
        
    def __iter__(self):
        return self
    
    def __next__(self):
        return self.next()

    def next(self):
        if self.isclosed: raise StopIteration()
        
        out = self.f.read(self.buff)
    
        if out == b'':
            self.close()
            raise StopIteration()
        else:
            if self.decode: return out.decode()
            else: return out
    
    def __del__(self):
        self.close()
    
    def close(self):
        self.isclosed = True
        self.f.close()