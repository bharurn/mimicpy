from pygmx.analysis import basic
from pygmx.run import gmxrun
import matplotlib.pyplot as plt
from tqdm import tqdm_notebook
import ipywidgets as widgets

def show(status=None):
    status = basic._checkstatus(status)
    
    temp = basic.Value.empty()
    pot = basic.Value.empty()

    for f in tqdm_notebook(status['edr'], desc='Loading Energies', leave=True):
        if f == 'em.edr': pass
        else:    
            pot += basic.energy(f, 'Potential')
            temp += basic.energy(f, 'Temperature')
        
    fig, ax = plt.subplots(2,2)
    #ax[0,0].set(xlabel=pot.xlabel, ylabel=pot.ylabel)
    ax[0,0].plot(pot.x, pot.y)
    ax[0,1].plot(temp.x, temp.y)
    plt.show()

def monitorLogs(log):
    
    df = basic.readLog(log)
    
    fig = plt.figure()
    ax = fig.add_subplot()
    line1, = ax.plot(df.reset_index()['Step'], df['Potential'])
    ax.set_autoscaley_on(True)
    plt.show()

    def plot_w(dataframe,xticker, yticker):
        line1.set_ydata(df[yticker])
        line1.set_xdata(df.reset_index()[xticker])
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw()
        fig.canvas.flush_events()

    widgets.interact(plot_w,
                     dataframe = widgets.fixed(df),
                 xticker = widgets.Dropdown(
            options=df.index.names,
            value='Step',
            description='X Axis:',
            disabled=False,
        ),
    yticker = widgets.Dropdown(
            options=df.columns,
            value='Potential',
            description='Y Axis:',
            disabled=False,
        )
    )

class GMXBoard(gmxrun.GMX):
    
    def __init__(self, gmx_run):
        self._status = gmx_run._status
    
    def show(self): show(self._status)
