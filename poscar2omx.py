#!/usr/bin/python3
import os
import numpy as np

if os.path.exists("mtd.md"):
    os.remove("mtd.md")

def write_file(file_name, content):
    with open(file_name, "w") as file:
        file.write(content)

def append_file(file_name, content):
    with open(file_name, "a") as file:
        file.write(content)

bohr2ang = 0.529177210903
ryd2har=0.5

nat=4
# nstep=1

cell,coords,energies=[],[],[]
with open("bilayer_graphene_mm/final.poscar", "r") as log_file:
    lines = log_file.readlines()
    for index, line in enumerate(lines):
        if (index>1) and (index<5):
            parts = line.split()
            vector = [float(parts[0]), float(parts[1]), float(parts[2])]
            cell.append(vector)
        elif (index>7) and (index<12):
            parts = line.split()
            vector = [float(parts[0]), float(parts[1]), float(parts[2])]
            coords.append(vector)

cell=np.array(cell)
coords=np.array(coords)
coords=np.dot(coords,cell)
print(cell.shape,coords.shape)

i=0
append_file("mtd.md",f"{nat}\n")
append_file("mtd.md",f"time=    {float(i):.4f} (fs) Energy= 0.000 (Hartree) Temperature=  300.000 (Given Temp.=  300.000) Cell_Vectors=  {cell[0][0]}  {cell[0][1]}  {cell[0][2]}  {cell[1][0]}  {cell[1][1]}  {cell[1][2]}  {cell[2][0]}  {cell[2][1]}  {cell[2][2]}\n")
append_file("mtd.md",f"C    {coords[0][0]}  {coords[0][1]}  {coords[0][2]}\n")
append_file("mtd.md",f"C    {coords[1][0]}  {coords[1][1]}  {coords[1][2]}\n")
append_file("mtd.md",f"C    {coords[2][0]}  {coords[2][1]}  {coords[2][2]}\n")
append_file("mtd.md",f"C    {coords[3][0]}  {coords[3][1]}  {coords[3][2]}\n")