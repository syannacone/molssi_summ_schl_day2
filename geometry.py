"""
geometry analysis functions
"""

import os
import numpy as np

def calculate_distance(rA, rB):
    """Calculate the distance between points A and B. Assumes rA and rB are numpy arrays."""
    dist_vec = (rA-rB)
    distance = np.linalg.norm(dist_vec)
    return distance

def calculate_distance_list(rA, rB):
    """Calculate the distance between points A and B. Assumes rA and rB are lists."""
    squared_sum = 0
    for dim in range(len(rA)):
        squared_sum += (rA[dim] - rB[dim])**2
    
    distance = np.sqrt(squared_sum)
    return distance


def build_bond_list(coordinates, max_bond=2.93, min_bond=0):
    """Builds list of bonds from atomic coordinates based on distance (short summary).

    Longer explanation of what the function does (extended summary).

    Parameters
    ----------
    coordinates : np.array
        An array of atomic coordinates. Size should be (n, 3), where n = # of particles
    max bond: float, optional parameter
        The max distance between atoms to be considered a bond. Default = 2.93 bohr
    min bond: float, optional parameter
        Minimum distance between atoms to be considered a bond

    Returns
    -------
    bonds : dict
        A dictionary of bonds with atom pair tuples as keys, and calculate bond lengths as values    
    """

    num_atoms = len(coordinates)
    
    bonds = {}
    
    for atom1 in range(num_atoms):
        for atom2 in range(atom1, num_atoms):
            distance = calculate_distance(coordinates[atom1], coordinates[atom2])
            
            if distance > min_bond and distance < max_bond:
                bonds[(atom1, atom2)] = distance 
    
    return bonds