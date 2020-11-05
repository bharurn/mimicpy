import re

def clean(txt, comments=None):
    if comments:  # Can be single string or a list of strings
        if isinstance(comments, str):
            comments = [comments]
        for c in comments:
            txt = re.sub(re.compile(str(c)+r"(.*)\n" ) ,"\n" , txt)  # Strip comments
    return re.sub(re.compile(r"[\n]+"), "\n", txt)

def print_dict(dct, col1, col2, printer):

    n1 = len(col1)+1
    n2 = len(col2)+1

    printer("+---------------------+")
    printer("| {:^{n1}}| {:^{n2}}|".format(col1, col2, n1=n1, n2=n2))
    printer("+---------------------+")
    for k,v in dct.items():
        printer("| {:^{n1}}| {:^{n2}}|".format(k, v, n1=n1, n2=n2))
        printer("+---------------------+")