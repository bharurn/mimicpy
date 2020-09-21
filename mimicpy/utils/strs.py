import inspect
import re


def f(s):
    frame = inspect.currentframe().f_back
    v = dict(**frame.f_globals)
    v.update(**frame.f_locals)
    return s.format(s, **v)


def clean(txt, comments=None):
    if comments:  # Comment can be a single string or a list of strings.
        if isinstance(comments, str):
            comments = [comments]

        for comment in comments:
            txt = re.sub(re.compile(f"{comment}(.*)\n" ) ,"\n" , txt)  # Strip comments.

    return re.sub(re.compile("\n{2,}" ) ,"\n" , txt)  # Removes double new lines.
