#!/usr/bin/python3
import os
import numpy as np

if os.path.exists("traj.md"):
    os.remove("traj.md")

def write_file(file_name, content):
    with open(file_name, "w") as file:
        file.write(content)

def append_file(file_name, content):
    with open(file_name, "a") as file:
        file.write(content)

bohr2ang = 0.529177210903
ryd2har=0.5

nat=4
nstep=14

cell,coords,energies=[],[],[]
with open("bilayer_graphene_mtd/log_qe.0.copy", "r") as log_file:
    lines = log_file.readlines()
    for index, line in enumerate(lines):
        if 'lattice vector [a.u.]' in line:
            vectors=[]
            for i in range(2, 5):
                vector_line = lines[index + i]
                parts = vector_line.split()
                if len(parts) == 4 and parts[0].isdigit():
                    vector = [float(parts[1]), float(parts[2]), float(parts[3])]
                    vectors.append(vector)
            cell.append(vectors)
        elif 'coord_xyz [a.u.]' in line:
            # print(index)
            # print(lines.index(line))
            coord=[]
            for i in range(2, nat+2):
                coord_line = lines[index + i]
                parts = coord_line.split()
                if len(parts) == 4 and parts[0].isdigit():
                    coord.append([float(parts[1]), float(parts[2]), float(parts[3])])
            coords.append(coord)
        elif '!    total energy' in line:
            parts=line.split()
            energies.append(float(parts[4]))
    cell=bohr2ang*np.array(cell)
    coords=bohr2ang*np.array(coords)
    # print(coords)
    energies=ryd2har*np.array(energies)
    # print(energies)
    # print(cell.shape,coords.shape,energies.shape)

for i in range(nstep):
    append_file("traj.md",f"{nat}\n")
    append_file("traj.md",f"time=    {float(i):.4f} (fs) Energy= {energies[i]} (Hartree) Temperature=  300.000 (Given Temp.=  300.000) Cell_Vectors=  {cell[i][0][0]}  {cell[i][0][1]}  {cell[i][0][2]}  {cell[i][1][0]}  {cell[i][1][1]}  {cell[i][1][2]}  {cell[i][2][0]}  {cell[i][2][1]}  {cell[i][2][2]}\n")
    append_file("traj.md",f"C    {coords[i][0][0]}  {coords[i][0][1]}  {coords[i][0][2]}\n")
    append_file("traj.md",f"C    {coords[i][1][0]}  {coords[i][1][1]}  {coords[i][1][2]}\n")
    append_file("traj.md",f"C    {coords[i][2][0]}  {coords[i][2][1]}  {coords[i][2][2]}\n")
    append_file("traj.md",f"C    {coords[i][3][0]}  {coords[i][3][1]}  {coords[i][3][2]}\n")