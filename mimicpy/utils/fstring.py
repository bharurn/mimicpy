import inspect

def f(s):
    frame = inspect.currentframe().f_back
    v = dict(**frame.f_globals)
    v.update(**frame.f_locals)
    return s.format(s, **v)