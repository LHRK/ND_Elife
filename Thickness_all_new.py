#!/usr/bin/env python
import argparse
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


parser = argparse.ArgumentParser(description='Measure minimum thickness for each lipid in both leaflets. The values are saved as two numpy arrays, which can begiven as input to the Thickness_apply_grid.py for distributing the values on a xy grid for plotning.')

#parser.add_argument('-f', dest='file_in' , action='store', type=str, help='-f for input xvg file')
parser.add_argument('-s', dest='gro', action='store', type=str,  help='-s for topology file')
parser.add_argument('-x ', dest='xtc', action='store', type=str,  help='-x for trajectory file')
parser.add_argument('--sel ', dest='selection', action='store', type=str,  help='--sel for selection for identifying the leaflet and calculating the thickness')
parser.add_argument('--cutoff ', dest='cut', action='store', default=15, type=float,  help='--cutoff for identification of the leaflets with MDAnalysis. Deafult is 15')

args = parser.parse_args()

gro    = args.gro
xtc    = args.xtc
sel    = args.selection
cutoff = args.cut

def get_thickness (sel1, sel2):
    loc_l = []
    for coor2 in sel2.positions[:,2]:
        loc_l.append(abs(sel1.positions[0][2]-coor2))
    l_idx = np.nonzero(np.array(loc_l))[0]
    l = np.array(loc_l)
    return np.min(l[l_idx]) #Return the minimum distance from sel1 to sel2, but not zero 

#### USER SPECIFIED #####
u = md.Universe(gro, xtc)
x_dim, y_dim, z_dim = u.dimensions[:3]
print ('Universe: {} {}'.format(gro,xtc))


nframes = len(u.trajectory[1:]) #Skipping the gro file
print ('Number of frames:', nframes)

L = LeafletFinder(u, sel, cutoff=cut)
leaflet0 = L.groups(0)
leaflet1 = L.groups(1)
top_leaf = " ".join(['{0:d}'.format(i) for i in leaflet0.resids ])
bot_leaf = " ".join(['{0:d}'.format(i) for i in leaflet1.resids ])


#For the top leaflet
grid_A = []

for f, ts in enumerate(u.trajectory[1:]):
    A = u.select_atoms('{1:s} and resid {0:s}'.format(top_leaf, sel), updating = True)
    B = u.select_atoms('{1:s} and resid {0:s}'.format(bot_leaf, sel), updating = True)
    nlipids = A.resids.shape[0]
    loc_array_thick  = np.zeros([nlipids])
    loc_array_x      = np.zeros([nlipids])
    loc_array_y      = np.zeros([nlipids])

    for nr, r in enumerate(A.resids): #Loop only over the top leaflet. 
        sel_loc = u.select_atoms('{0:s} and resid {1:d}'.format(sel, r))
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
    A = u.select_atoms('{1:s} and resid {0:s}'.format(top_leaf, sel), updating = True)
    B = u.select_atoms('{1:s} and resid {0:s}'.format(bot_leaf, sel), updating = True)
    nlipids = B.resids.shape[0]
    loc_array_thick  = np.zeros([nlipids])
    loc_array_x      = np.zeros([nlipids])
    loc_array_y      = np.zeros([nlipids])

    for nr, r in enumerate(B.resids): #Loop only over the top leaflet. 
        sel_loc = u.select_atoms('{0:s} and resid {1:d}'.format(sel, r))
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

