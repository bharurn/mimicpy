def read(file, mode):
    with open(file, mode) as f:
        out = f.read()
    return out

def write(out, file, mode):
    with open(file, mode) as f:
        f.write(out)