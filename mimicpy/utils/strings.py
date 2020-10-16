import re

def clean(txt, comments=None):
    if comments:  # Can be single string or a list of strings
        if isinstance(comments, str):
            comments = [comments]
        for c in comments:
            txt = re.sub(re.compile(f"{c}(.*)\n" ) ,"\n" , txt)  # Strip comments
    return re.sub(re.compile(r"[\n]+"), "\n", txt)