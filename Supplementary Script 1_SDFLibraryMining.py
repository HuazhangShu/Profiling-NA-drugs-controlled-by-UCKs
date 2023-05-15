"""
Author: Huazhang Shu
Date: Jan., 19th, 2023

This module reads .sdf files and split the library into each molecule with four attributes:  (1) name, 
(2) formula, (3) CAS, (4) coordinates. Further processing was done by checking the intersection between two 
libraries. Note that the library itself may lack some information of some molecules. 

The intersection was defined by searching for common CAS numbers existing in both two libraries.
"""

# Library preparation
import os


# Class initialization
class Molecule:
    def __init__(self, name='nan', formula='nan', CAS='nan', coordinates=[]):
        self.name = name
        self.formula = formula
        self.CAS = CAS
        self.coordinates = coordinates


# The following functions read split text containing information of one molecule and return corresponding information
# as names suggest.
def get_name(file):
    with open(file, 'r') as f0:
        if_name = False
        mol_name = 'nan'
        for f0_line in f0:
            if '<Name>' in f0_line:
                if_name = True
            elif if_name:
                mol_name = f0_line.strip('\n')
                if_name = False
    return mol_name


def get_formula(file):
    with open(file, 'r') as f0:
        if_formula = False
        mol_formula = 'nan'
        for f0_line in f0:
            if '<Formula>' in f0_line:
                if_formula = True
            elif if_formula:
                mol_formula = f0_line.strip('\n')
                if_formula = False
    return mol_formula


def get_CAS(file):
    with open(file, 'r') as f0:
        if_CAS = False
        mol_CAS = 'nan'
        for f0_line in f0:
            if '<CAS>' in f0_line:
                if_CAS = True
            elif if_CAS:
                mol_CAS = f0_line.strip('\n')
                if_CAS = False
    return mol_CAS


def get_coordinates(file):
    start = 0
    end = 0
    with open(file, 'r') as f0:
        f0_lines = f0.readlines()
        mol_csChFnd80 = 'nan'
        for i, f0_line in enumerate(f0_lines):
            if 'csChFnd80' in f0_line:
                start = i
            if 'END' in f0_line:
                end = i
            mol_csChFnd80 = f0_lines[start: end + 1]
    return mol_csChFnd80


# This function split .sdf library and extract each molecule as an object stored in a molecule list.
# Molecules order will be preserved in the same way as in the .sdf file.
def get_molecules(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        mol_list = []
        last_end_index = 0
        for i, line in enumerate(lines):
            if "$$$$" in line:
                end_index = i
                new_lines = lines[last_end_index: end_index]
                with open('output.txt', 'w') as f2:
                    for new_line in new_lines:
                        f2.write(new_line)
                mol_list.append(Molecule(name=get_name('output.txt'),
                                         formula=get_formula('output.txt'),
                                         CAS=get_CAS('output.txt'),
                                         coordinates=get_coordinates('output.txt')))
                last_end_index = end_index + 1
    os.remove('output.txt')
    return mol_list


if __name__ == "__main__":
    # Create molecule objects for each entry in the nucleoside analogue library and FDA-approved drug library
    # Write each CAS number into different files to store
    na_list = get_molecules('20230106-L7200-Nucleoside-Analogue-Library.SDF') # Change the file name for new dataset
    with open('na_cas_list.txt', 'w') as f1:
        for mol in na_list:
            f1.write(mol.CAS + '\n')
    fda_list = get_molecules('20221216-L1300-FDA-approved-Drug-Library.SDF') # Change the file name for new dataset
    with open('fda_cas_list.txt', 'w') as f2:
        for mol in fda_list:
            f2.write(mol.CAS + '\n')

    # See intersections
    with open('na_cas_list.txt', 'r') as f1:
        with open('fda_cas_list.txt', 'r') as f2:
            list_na = [i.strip('\n') for i in f1.readlines()]
            list_fda = [i.strip('\n') for i in f2.readlines()]

    with open('intersection.txt', 'w') as f3:
        for i in list_na:
            if i in list_fda:
                if i == 'nan':
                    continue
                f3.write(i + '\n')

    na_CAS_list = [i.CAS for i in na_list]
    na_coord_list = [i.coordinates for i in na_list]
    na_name_list = [i.name for i in na_list]

    # Derive a complete table with molecule names and CAS numbers
    with open('intersection.txt', 'r') as f1:
        intersection_list = f1.readlines()

    with open('intersection_full.txt', 'w') as f2:
        for mol in intersection_list:
            index = na_CAS_list.index(mol.strip('\n'))
            f2.write(na_CAS_list[index] + '\t' + na_name_list[index] + '\n')
