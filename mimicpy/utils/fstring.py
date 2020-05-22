# This is part of MiMiCPy

"""

This module contains the f function to simulate
f-string behavior in python version 3.6-

"""

import inspect

def f(s):
    frame = inspect.currentframe().f_back
    v = dict(**frame.f_globals)
    v.update(**frame.f_locals)
    return s.format(s, **v)