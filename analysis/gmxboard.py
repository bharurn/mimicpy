from pygmx.analysis import basic
from pygmx.run import gmxrun
from pygmx import host
import matplotlib.pyplot as plt


class GMXBoard(gmxrun.GMX):
    
    def show():
        file = host.cmd.ls(file_eval=lambda a: True if a.startswith('md') and a.endswith('.edr') else False)
        
        fig, axs = plt.subplots(2, 2)
        axs[0, 0].plot(x, y)
        axs[0, 0].set_title('Axis [0, 0]')
        axs[0, 1].plot(x, y, 'tab:orange')
        axs[0, 1].set_title('Axis [0, 1]')
        axs[1, 0].plot(x, -y, 'tab:green')
        axs[1, 0].set_title('Axis [1, 0]')
        axs[1, 1].plot(x, -y, 'tab:red')
        axs[1, 1].set_title('Axis [1, 1]')

        for ax in axs.flat:
            ax.set(xlabel='x-label', ylabel='y-label')
            