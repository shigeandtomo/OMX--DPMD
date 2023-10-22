#!/usr/bin/python3

import os,sys

os.chdir("10.data/")
# 2回目以降
# os.remove("met.evp")
# os.remove("met.for")
# os.remove("met.in")
# os.remove("met.pos")

Ang2Bohr=1.88973 #Ang->Bohr
force=1.0 #Har/Bohr->Har/Bohr

nat=5
ntyp=2

# ファイルの書き込み関数
def write_file(file_name, content):
    with open(file_name, 'w') as file:
        file.write(content)

# ファイルの追加書き込み関数
def append_file(file_name, content):
    with open(file_name, 'a') as file:
        file.write(content)

# //iterout.c
# /* 1: */
# /* 2,3,4: */
# /* 5,6,7: force *   
# /* 8: x-component of velocity */
# /* 9: y-component of velocity */
# /* 10: z-component of velocity */
# /* 11: Net charge, electron charge is defined to be negative. */
# /* 12: magnetic moment (muB) */
# /* 13,14: angles of spin */  

# met.evp ファイルの生成
write_file("met.evp", "#   nfi    time(ps)        ekinc        T_cell(K)     Tion(K)             etot\n")
with open("met.md", 'r') as md_file:
    lines = md_file.readlines()
    nfi=0
    for index, line in enumerate(lines):
        if "time" in line:
            nfi+=1
            parts = line.split()
            time_ps = float(parts[1]) * 1e-3 #fs->ps
            ekinc = 0.0
            T_cell = 0.0
            Tion = float(parts[7])
            etot = float(parts[4]) #Hartree->Hartree
            evp_line = f"     {nfi}  {time_ps:.6e}  {ekinc:.6e}  {T_cell:.6e}  {Tion:.6e}        {etot:.6e}\n"
            append_file("met.evp", evp_line)

# met.for ファイルの生成
with open("met.md", 'r') as md_file:
    lines = md_file.readlines()
    cnt=0
    for index, line in enumerate(lines):
        if "C " in line or "H " in line:
            cnt+=1
            if (cnt -1)% nat == 0:
                for_line = f"{      cnt // nat + 1}\n"
                append_file("met.for",for_line)
            parts = line.split()
            for_line = f"    {float(parts[4])*force:.6e}    {float(parts[5])*force:.6e}    {float(parts[6])*force:.6e}\n"
            append_file("met.for", for_line)

# met.pos ファイルの生成
with open("met.md", 'r') as md_file:
    lines = md_file.readlines()
    cnt=0
    for index, line in enumerate(lines):
        if "C " in line or "H " in line:
            cnt+=1
            if (cnt -1)% nat == 0:
                for_line = f"{cnt // nat + 1}\n"
                append_file("met.pos",for_line)
            parts = line.split()
            for_line = f"    {float(parts[1])*Ang2Bohr:.6e}    {float(parts[2])*Ang2Bohr:.6e}    {float(parts[3])*Ang2Bohr:.6e}\n"
            append_file("met.pos", for_line)

# met.in ファイルの生成
write_file("met.in", f" &control\n /\n &system\n    ibrav = 14,\n    celldm(1) = 12.0,\n    celldm(2) = 1.0,\n    celldm(3) = 1.0,\n    celldm(4) = 0.0,\n    celldm(5) = 0.0,\n    celldm(6) = 0.0,\n    nat  = {nat},\n    ntyp = {ntyp},\n /\n &electrons\n /\n &ions\n /\n &cell\n /\nATOMIC_SPECIES\n")

with open("Methane.dat", 'r') as dat_file:
    lines = dat_file.readlines()
    species = []
    species_mode=False
    for line in lines:
        if "<Definition.of.Atomic.Species" in line:
            species_mode = True
        elif "Definition.of.Atomic.Species>" in line:
            species_mode = False
        elif species_mode:
            parts = line.split()
            species.append(f"{parts[0]} 1.00d0 H.blyp-vbc.UPF")
    species_content = "\n".join(species)
    append_file("met.in", species_content)
    append_file("met.in","\n")

with open("Methane.dat", 'r') as dat_file:
    lines = dat_file.readlines()
    coords = []
    coords_mode=False
    for line in lines:
        if "<Atoms.SpeciesAndCoordinates" in line:
            coords_mode = True
        elif "Atoms.SpeciesAndCoordinates>" in line:
            coords_mode = False
        elif coords_mode:
            parts = line.split()
            coords.append(f"{parts[1]} {float(parts[2])*Ang2Bohr} {float(parts[3])*Ang2Bohr} {float(parts[4])*Ang2Bohr}")
    coords_content = "\n".join(coords)
    append_file("met.in", "ATOMIC_POSITIONS (bohr)\n" + coords_content)
