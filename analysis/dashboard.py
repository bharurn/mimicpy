import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display

def show(*args, **font):
    fnt = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 20}
    
    if font is not None: fnt.update(font)
    
    plt.rc('font', **fnt)

    tab = widgets.Tab()
    chld = []
    
    for i,a in enumerate(args):
        tab.set_title(i, a.title)
        chld.append(a.getBox())
    
    tab.children = chld
    display(tab)