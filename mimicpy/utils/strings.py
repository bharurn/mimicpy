import re

def clean(string, comments=None):
    if comments:  # comments can be a single string or a list of strings.
        if isinstance(comments, str):
            comments = [comments]
        for comment in comments:
            string = re.sub(re.compile(f"{comment}(.*)\n" ) ,"\n" , string)  # Strip comments.
    return re.sub(re.compile("\n{2,}" ) ,"\n" , string)  # Remove double new lines.
