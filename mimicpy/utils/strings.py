import re

def clean(txt, comments=None):
    if comments:  # Can be single string or a list of strings
        if isinstance(comments, str):
            comments = [comments]
        for c in comments:
            txt = re.sub(re.compile(str(c)+r"(.*)\n" ) ,"\n" , txt)  # Strip comments
    return re.sub(re.compile(r"[\n]+"), "\n", txt)

def print_table(dct, printer):
    cols = dct.keys()
    n = [len(c)+1 for c in cols]
    template = "| " + "| ".join(["{:^" + str(i) + "}" for i in n]) + "|"
    
    dashes = '-'*(sum(n)+2*len(n)-1)
    
    printer("+{}+".format(dashes))
    printer(template.format(*cols))
    printer("+{}+".format(dashes))
    
    vals = dct.values()
    lst = list(map(list, zip(*vals))) # transpose list
    
    for i in lst:
        printer(template.format(*i))
        printer("+{}+".format(dashes))

def print_dict(dct, col1, col2, printer):
    new_dct = {col1: list(dct.keys()), col2: list(dct.values())}
    print_table(new_dct, printer)