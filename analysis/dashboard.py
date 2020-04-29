from itertools import cycle
import random

try:
    import matplotlib.pyplot as plt
    import ipywidgets as widgets
    from IPython.display import display
    from ipywidgets import interactive, Layout, HBox, VBox
except ImportError:
    raise Exception("Please install the analysis verions of the package")
            

#decorator syntax:
#@plotBoxGen(title of tab, args to be passed to func..)
def plotBoxGen(*args):
    def decorator(function):
        def wrapper(flag, *args1):
            if flag:
                a = list(args)
                return (a[0], a[1:])
            else: return function(*args1)
        return wrapper
    return decorator

class PlotBox:
    def __init__(self, df_gen, *args):
        self.line1 = None
        self.ax = None
        self.fig = None
        self.df_gen = df_gen
        self.df_args = args
        self.df = None
        self.df_no_ndx = None
        self.title, self.index = self.df_gen(True)
        self.xdef_opt = "--Select to Calculate--"
        self.ydef_opt = "--Not Calculated--"
       
    def plot(self, x, y):
        if self.df is None:
            self.df = self.df_gen(False, *self.df_args)
            self.df_no_ndx = self.df.reset_index()
            self.fixColor()
            self.fig = plt.figure()
            
        if y == self.ydef_opt:
            self.y.options = self.df.columns
            val = x
            self.x.options = self.index
            self.x.value = val
            return
            
        if x == self.xdef_opt or y == self.ydef_opt: return
        
        else:   
            plt.close()
            fig = plt.figure(num=f"{y} vs {x} | Use Controls Above to Change Axis & Below to Interact", figsize=(4.5,3))
            ax = fig.add_subplot(111)
            ax.plot(self.df_no_ndx[x], self.df[y], color=self.colors[self.df.columns.get_loc(y)],\
                             marker='o', linewidth=2.5, markersize=8, linestyle='-.')
            ax.set_autoscaley_on(True)
            ax.set_ylabel(y)
            ax.set_xlabel(x)
            plt.show()
       
    def fixColor(self):
        c_list = ["r","b","g","c","m","y","k"]
        if len(c_list) > len(self.df.columns):
            self.colors = random.sample(["r","b","g","c","m","y","k"], len(self.df.columns))
        else:
            self.colors = []
            for i,a in enumerate(cycle(c_list)):
                if i >= len(self.df.columns): break
                else: self.colors.append(a)
        
    def getDropDown(self, labels, val, desc):
        labels.append(val)
        
        return widgets.Dropdown(
            options=labels,
            value=val,
            description=desc
        )
        
    def getBox(self):
        
        self.box = VBox(children = self.getWidget().children)
        self.box.layout = Layout(margin='23px 0px 0px 100px')
        self.box.children[0].layout = Layout(left='0px', top="0px")
        self.box.children[1].layout = Layout(left='350px', top="-32px")
        
        return self.box
    
    def getHBox(self):
        
        self.box = HBox(children = self.getWidget().children)
        
        return self.box
    
    def getWidget(self):
        self.x = self.getDropDown(self.index.copy(), self.xdef_opt, 'X Axis:')
        self.y = self.getDropDown([], self.ydef_opt, 'Y Axis:')
        return interactive(lambda x, y: self.plot(x, y), x = self.x, y = self.y)
    
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