from itertools import cycle
import random
import threading
import time
import sys
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display
from ipywidgets import interactive, Layout, HBox, VBox

#decorator PlotBoxDF
class _PlotBoxDF(object):
    def __init__(self, foo, x_axis=['Step', 'Time']):
        self.foo = foo
        self.x_axis = x_axis
        self.estimated_time = 2 # not loading bar
        self.tstep = 0.2
    
    def pbar_wrapper(self, function, args=[], kwargs={}):
        ret = [None]  # Mutable var so the function can store its return value
        def myrunner(function, ret, *args, **kwargs):
            ret[0] = function(*args, **kwargs)

        thread = threading.Thread(target=myrunner, args=(function, ret) + tuple(args), kwargs=kwargs)
        pbar = tqdm(total=self.estimated_time, file=sys.stdout,\
                bar_format='{desc}:{bar}', desc='Calculating', leave=False, disable=False)

        thread.start()
        while thread.is_alive():
            thread.join(timeout=self.tstep)
            pbar.update(self.tstep)
    
        pbar.close()
    
        return ret[0]

    def __call__(self, *args, **kwargs):
        if self.estimated_time != 0:
            return self.pbar_wrapper(self.foo, args, kwargs)
        else:
            print(args)
            return self.foo(*args, **kwargs)
    
# wrap _PlotBoxDF to allow for deferred calling
def PlotBoxDF(function=None, x_axis=['Time', 'Step']):
    if function:
        return _PlotBoxDF(function)
    else:
        def wrapper(function):
            return _PlotBoxDF(function, x_axis=x_axis)

        return wrapper
            
class StaticPlot:
    def __init__(self, df_gen, *args):
        self._df_gen = df_gen
        self._df_args = args
        self._df = None
        self._df_no_ndx = None
        self.title = 'Plot'
        self._index = self._df_gen.x_axis
        self._xdef_opt = "--Select to Calculate--"
        self._ydef_opt = "--Not Calculated--"
        self.figsize = (7,5.5)
        self.ax_kwargs = {'marker':'o', 'linewidth':2, 'markersize':4.5, 'linestyle':'-.'}
    
    def setAxisProps(self, **kwargs): self.ax_kwargs.update(kwargs)
    
    def _init(self):
        self._df = self._df_gen(*self._df_args)
        self._df_no_ndx = self._df.reset_index()
        
    def _plot(self, x, y):
        if self._df is None:
            self._init()
            self._fixColor()
            
        if y == self._ydef_opt:
            self._y.options = self._df.columns
            val = x
            self._x.options = self._index
            self._x.value = val
            return
            
        if x == self._xdef_opt or y == self._ydef_opt: return
        
        else:   
            plt.close()
            fig = plt.figure(num=f"{y} vs {x}", figsize=self.figsize)
            ax = fig.add_subplot(111)
            ax.plot(self._df_no_ndx[x], self._df[y], color=self._colors[self._df.columns.get_loc(y)], **self.ax_kwargs)
            ax.set_autoscaley_on(True)
            ax.set_ylabel(y)
            ax.set_xlabel(x)
            fig.tight_layout()
            plt.show()
       
    def _fixColor(self):
        c_list = ["r","b","g","c","m","y","k"]
        if len(c_list) > len(self._df.columns):
            self._colors = random.sample(["r","b","g","c","m","y","k"], len(self._df.columns))
        else:
            self._colors = []
            for i,a in enumerate(cycle(c_list)):
                if i >= len(self._df.columns): break
                else: self._colors.append(a)
        
    def _getDropDown(self, labels, val, desc):
        labels.append(val)
        
        return widgets.Dropdown(
            options=labels,
            value=val,
            description=desc
        )
        
    def getBox(self):
        
        self._box = VBox(children = self.getWidget().children)
        self._box.layout = Layout(margin='23px 0px 0px 40px')
        self._box.children[0].layout = Layout(left='20px', top="0px")
        self._box.children[1].layout = Layout(left='450px', top="-32px")
        
        return self._box
    
    def getHBox(self):
        
        self._box = HBox(children = self.getWidget().children)
        
        return self._box
    
    def getWidget(self):
        self._x = self._getDropDown(self._index.copy(), self._xdef_opt, 'X Axis:')
        self._y = self._getDropDown([], self._ydef_opt, 'Y Axis:')
        return interactive(lambda x, y: self._plot(x, y), x = self._x, y = self._y)

class MonitorPlot(StaticPlot):
    def __init__(self, df_gen, update, *args):
        super().__init__(df_gen, *args)
        self.update = update
        self._time = time.time()
    
    def _plot(self, x, y):
        if (time.time() - self._time) >= self.update:
            self._init()
            self._time = time.time()
    
        super()._plot(x, y)
        
def show(*args, **font):
    fnt = {'family' : 'Arial',
        'size'   : 14}
    
    if font is not None: fnt.update(font)
    
    plt.rc('font', **fnt)

    tab = widgets.Tab()
    chld = []
    
    for i,a in enumerate(args):
        tab.set_title(i, a.title)
        chld.append(a.getBox())
    
    tab.children = chld
    display(tab)
    
def watch(gmx_log=None, cpmd_log=None):
    parse = sys.modules[__package__ + '.parse']
    if gmx_log:
        p1 = MonitorPlot(parse.log, 30, gmx_log)
        p1.title = 'Gromacs Log'
        show(p1)