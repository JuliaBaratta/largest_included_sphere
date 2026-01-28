import numpy as np
np.set_printoptions(suppress=True, precision=4)
from pathlib import Path
from ase import io
from ase.visualize import view
from scipy.spatial import cKDTree
from ase import Atom
import json

def lonely_point(path_to_file):
    '''
    Function that finds the lonliest grid point in every crystal and puts a atom there
    
    Inputs: 
        - path_to_file (string): path to the .cif file of the crystal

    Outputs: 
        - .cif file with the largest included sphere added

    '''
    #load in file
    cif = io.read(path_to_file)
    atom_positions = cif.get_positions()
    cell = cif.get_cell()

    #tree for fast sorting
    tree = cKDTree(atom_positions)

    # Create nxnxn grid to search structure by 
    n_points = 100
    frac_x = np.linspace(0, 1, n_points)
    frac_y = np.linspace(0, 1, n_points)
    frac_z = np.linspace(0, 1, n_points)

    #puts the nxnxn grid onto a mesh and indexes it for faster sorting
    frac_grid = np.meshgrid(frac_x, frac_y, frac_z, indexing='ij')
    fractional_coords = np.column_stack([grid.ravel() for grid in frac_grid])
    grid_cartesian = fractional_coords @ cell.array

    d_matrix, nearest_atom_indices = tree.query(grid_cartesian)

    # Sort
    sorted_indices = np.argsort(d_matrix)[::-1]
    loneliest_idx = sorted_indices[0]

    # Option to print out the sorted indicies for more control over where you put the 'X' atom
    # with open("mydata.json", "w") as file:
	  #    json.dump(fractional_coords[sorted_indices].tolist(), file)

    #some potentially helpful info to print out about your structure
    print(f"\n=== LONELIEST POINT ===")
    print(f"Distance to nearest atom: {d_matrix[loneliest_idx]:.4f} Å")
    print(f"Fractional position: {fractional_coords[loneliest_idx]}")
    print(f"Cartesian position: {grid_cartesian[loneliest_idx]} Å")
    print(f"Nearest atom index: {nearest_atom_indices[loneliest_idx]}")
    print(f"Nearest atom type: {cif[nearest_atom_indices[loneliest_idx]].symbol}")
    
    # copy structure and add lonely point
    dummy_structure = cif.copy()
    loneliest_point_cartesian_coordinates = fractional_coords[loneliest_idx] @ cif.get_cell().array
  
    # option to add a second lonliest point with a specific index
    #second_loneliest_point_cartesian_coordinates = [0.494949494949495, 0.494949494949495, 0.494949494949495] @ cif.get_cell().array

    # Append the dummy atom
    dummy_structure.append(Atom('X', loneliest_point_cartesian_coordinates))

    # option to add that second point
    #dummy_structure.append(Atom('X', second_loneliest_point_cartesian_coordinates))

    view(dummy_structure)
    # commenting out for rn, don't want to view every structure in the mp
    io.write(f'lis_{path_to_file}', dummy_structure)

#replace file.cif with the name of your file
lonely_point(Path('./file.cif'))
