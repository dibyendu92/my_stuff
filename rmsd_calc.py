#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

###PDB Record Format###

# COLUMNS        DATA TYPE     FIELD        DEFINITION
# -------------------------------------------------------------------------------------
# 1 - 6        Record name   "ATOM  "
# 7 - 11       Integer       serial       Atom serial number.
# 13 - 16      Atom          name         Atom name.
# 17           Character     altLoc       Alternate location indicator.
# 18 - 20      Residue name  resName      Residue name.
# 22           Character     chainID      Chain identifier.
# 23 - 26      Integer       resSeq       Residue sequence number.
# 27           AChar         iCode        Code for insertion of residues.
# 31 - 38      Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
# 39 - 46      Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
# 47 - 54      Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
# --------------------------------------------------------------------------------------

def calc_COM(atom_dic):
    mass_dic = {"N":14.0, "O": 16.0, "S": 32.0, "H": 1.0, "C": 12.0}
    COM_x = 0.0
    COM_y = 0.0
    COM_z = 0.0
    mass_total = 0.0
    ## atom_dic format
    ## atom[Atom serial number] = [Atom name,Residue name, Residue seq number, Coord X, Coord Y, Coord Z] 
    for keys in atom_dic:
        if atom_dic[keys][0].strip()[0] in mass_dic:
           mass_total = mass_total + mass_dic[atom_dic[keys][0].strip()[0]] 
           COM_x = COM_x + mass_dic[atom_dic[keys][0].strip()[0]]*atom_dic[keys][3]
           COM_y = COM_y + mass_dic[atom_dic[keys][0].strip()[0]]*atom_dic[keys][4]
           COM_z = COM_z + mass_dic[atom_dic[keys][0].strip()[0]]*atom_dic[keys][5]
        else:
           print "unknown atom name", atom_dic[keys][0].strip()[0]
        
    return (COM_x/mass_total, COM_y/mass_total, COM_z/mass_total)

def find_dist(pdb_file1, pdb_file2):
    atom1 = {}
    atom2 = {}
    res_name = ["GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "LYS", "HIS", "PHE", "TYR", "TRP", "ASN", "GLN", "SER", "THR", "ARG", "ASP", "GLU", "CYS", "MET", "CYI", "HIP"]
    with open(pdb_file1, "r") as input_file1:
         for line in input_file1:
             if line[0:6] == "ATOM  ":
                atom1[line[6:11]] = [line[12:16],line[17:20], line[22:26], float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
#               atom[Atom serial number] = [Atom name,Residue name, Residue seq number, Coord X, Coord Y, Coord Z] 
    (COM_X1, COM_Y1, COM_Z1) = calc_COM(atom1) 
    with open(pdb_file2, "r") as input_file2:
         for line in input_file2:
             if line[0:6] == "ATOM  ":
                atom2[line[6:11]] = [line[12:16],line[17:20], line[22:26], float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
#               atom[Atom serial number] = [Atom name,Residue name, Residue seq number, Coord X, Coord Y, Coord Z]
    (COM_X2, COM_Y2, COM_Z2) = calc_COM(atom2) 
    if len(atom1) == len(atom2):
         res_dist_list = {}
         total_dist2 = 0.0
         for keys in atom1:
            if atom2[keys][0]+atom2[keys][1]+atom2[keys][2] == atom1[keys][0]+atom1[keys][1]+atom1[keys][2]:
               if atom1[keys][1] in res_name: 
                  a = np.array((atom1[keys][3], atom1[keys][4], atom1[keys][5]))
                  b = np.array((atom2[keys][3], atom2[keys][4], atom2[keys][5])) 
                  dist = np.linalg.norm(a-b)
                  dist2 = dist*dist
                  total_dist2 +=dist2
                  if atom1[keys][2].strip() not in res_dist_list:
                     res_dist_list[atom1[keys][2].strip()] = [dist*dist]
                  else:
                     res_dist_list[atom1[keys][2].strip()].append(dist*dist)
         rmsd_res = {}
         for key in res_dist_list:
             rmsd_res[key] = 0.0
             for i in range(len(res_dist_list[key])):
                 rmsd_res[key] +=res_dist_list[key][i]
             rmsd_res[key] =np.sqrt(rmsd_res[key]/(len(res_dist_list[key])-1))
     
         if sys.argv[3] != "FULL":
             x=[]
             y=[]
             for key in rmsd_res:
                 x.append(key)
                 y.append(rmsd_res[key])
             plt.bar(x,y)
             plt.show()
         else:     
             print(np.sqrt(total_dist2/len(atom1)))     
    return None
          
pdb_file1 = sys.argv[1]
pdb_file2 = sys.argv[2]

find_dist(pdb_file1, pdb_file2)
