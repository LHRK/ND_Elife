#!/usr/bin/env python

import numpy as np
import MDAnalysis as md
import matplotlib.pyplot as plt 
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from numpy import ma
from matplotlib import ticker, cm


#######################################################
# Loads in the files from the script Thickess_all.py  #
# Divides the values unto a grid                      #
#######################################################


grid_t_A    = np.load('Thickness_Top_all.npy') #Has the shape [nframes, nlipids, 3] 
grid_t_B    = np.load('Thickness_Bot_all.npy')

u = md.Universe('topology.pdb', 'all_simulations.xtc')
x_dim, y_dim, z_dim = u.dimensions[:3]

nbins = 5

x_steps = np.linspace(0, x_dim, num=nbins+1)
y_steps = np.linspace(0, y_dim, num=nbins+1)

nframes = len(u.trajectory[1:]) #Skip gro file

def apply_grid(grid_t, nframes, nbins):
    grid = np.zeros([nbins,nbins, nframes])
    
    for f in range(nframes):
        for i in range(nbins):
            for j in range(nbins):
                x = grid_t[f,:,0]
                y = grid_t[f,:,1]
    
                bool_idx = [(x>=x_steps[i]) & (x<x_steps[i+1]) & (y>=y_steps[j]) & (y<y_steps[j+1])]
                idx = np.where(bool_idx[0] ==True)
                ava_bin = np.average(grid_t[f,idx,2])
                
                if np.isnan(ava_bin):
                    grid[i,j,f] = 0
                else:
                    grid[i,j,f] = ava_bin
    return grid

grid_A = apply_grid(grid_t_A, nframes, nbins)
grid_B = apply_grid(grid_t_B, nframes, nbins)

grid_ava = (grid_A + grid_B) / 2

np.save('Thickness_Top_grid_{0:d}.npy'.format(nbins), grid_A)
np.save('Thickness_Bot_grid_{0:d}bins.npy'.format(nbins), grid_B)
np.save('Thickness_averaged_grid_{0:d}bins.npy'.format(nbins), grid_ava)
print ('Files saved as: Thickness_Top_grid_{0:d}.npy, Thickness_Bot_grid_{0:d}.npy, Thickness_averaged_grid_{0:d}.npy'.format(nbins))


