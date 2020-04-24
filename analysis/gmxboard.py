from pygmx.analysis import basic
from pygmx.run import gmxrun
import matplotlib.pyplot as plt


class GMXBoard(gmxrun.GMX):
    
    def show(self):
        temp = basic.Value.empty()
        pot = basic.Value.empty()

        for f in self.gethistory('edr'):
            if f == 'em.edr':
                pass
            else:    
                print(f"Reading {f}")
                pot += basic.energy(f, 'Potential')
                temp += basic.energy(f, 'Temperature')
        
        fig, ax = plt.subplots(2,2)
        #ax[0,0].set(xlabel=pot.xlabel, ylabel=pot.ylabel)
        ax[0,0].plot(pot.x, pot.y)
        ax[0,1].plot(temp.x, temp.y)
        plt.show()