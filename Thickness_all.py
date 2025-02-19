#!/usr/bin/env python

import numpy as np
import MDAnalysis as md
from MDAnalysis.analysis.leaflet import LeafletFinder
from numpy import ma

##########################################################
# Calculates the minimum thickness for each              #
# lipid in the top and bottom leaflet                    #
# The thickness values and corresponding x,y coordinates #
# are then saved to files                                #
##########################################################


def get_thickness (sel1, sel2):
    loc_l = []
    for coor2 in sel2.positions[:,2]:
        loc_l.append(abs(sel1.positions[0][2]-coor2))
    l_idx = np.nonzero(np.array(loc_l))[0]
    l = np.array(loc_l)
    return np.min(l[l_idx]) #Return the minimum distance from sel1 to sel2, but not zero 

#### USER SPECIFIED #####
u = md.Universe('topology.pdb', 'all_simulations.xtc')
x_dim, y_dim, z_dim = u.dimensions[:3]
print ('Universe loaded')


nframes = len(u.trajectory[1:]) #Skipping the gro file
print ('Number of frames:', nframes)

L = LeafletFinder(u, 'resname DMPC and name P')
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)
top_leaf = " ".join(['{0:d}'.format(i) for i in leaflet0.resids ])
bot_leaf = " ".join(['{0:d}'.format(i) for i in leaflet1.resids ])


#For the top leaflet
grid_A = []

for f, ts in enumerate(u.trajectory[1:]):
    A = u.select_atoms('resid {0:s} and name P'.format(top_leaf), updating = True)
    B = u.select_atoms('resid {0:s} and name P'.format(bot_leaf), updating = True)
    nlipids = A.resids.shape[0]
    loc_array_thick  = np.zeros([nlipids])
    loc_array_x      = np.zeros([nlipids])
    loc_array_y      = np.zeros([nlipids])

    for nr, r in enumerate(A.resids): #Loop only over the top leaflet. 
        sel_loc = u.select_atoms('resname DMPC and name P and resid {0:d}'.format(r))
        thick = get_thickness(sel_loc, B)
        x,y,z = sel_loc.positions[0]
        loc_array_x[nr] = x
        loc_array_y[nr] = y
        loc_array_thick[nr] = thick
    grid_loc = np.column_stack((loc_array_x, loc_array_y, loc_array_thick))
    grid_A.append(grid_loc)

grid_t_A = np.array(grid_A)

# For the bottom leaflet
grid_B = []

for f, ts in enumerate(u.trajectory[1:]):
    A = u.select_atoms('resid {0:s} and name P'.format(top_leaf), updating = True)
    B = u.select_atoms('resid {0:s} and name P'.format(bot_leaf), updating = True)
    nlipids = B.resids.shape[0]
    loc_array_thick  = np.zeros([nlipids])
    loc_array_x      = np.zeros([nlipids])
    loc_array_y      = np.zeros([nlipids])

    for nr, r in enumerate(B.resids): #Loop only over the top leaflet. 
        sel_loc = u.select_atoms('resname DMPC and name P and resid {0:d}'.format(r))
        thick = get_thickness(sel_loc, A)
        x,y,z = sel_loc.positions[0]
        loc_array_x[nr] = x
        loc_array_y[nr] = y
        loc_array_thick[nr] = thick
    grid_loc = np.column_stack((loc_array_x, loc_array_y, loc_array_thick))
    grid_B.append(grid_loc)

grid_t_B = np.array(grid_B)

np.save('Thickness_Top_all.npy', grid_t_A)
np.save('Thickness_Bot_all.npy', grid_t_B)
print ('Files saved as: Thickness_Top_all.npy and Thickness_Bot_all.npy for the top and bottom leaflet')

