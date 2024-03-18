#!/bin/bash

cd OMX
cat >h2o.evp<<EOF
#   nfi    time(ps)        ekinc        T_cell(K)     Tion(K)             etot            
EOF
grep 'time' H2O.md |awk '{printf "%d %.6e %.6e %.6e %.6e %.6e\n", NR, $2*1e-3, "0", "0",  $8, $5}' >>h2o.evp
grep 'O\|H ' H2O.md |awk '(NR-1)%3==0 {printf "%.0f\n", NR/3+1} {print $8*1.99467e-2, $9*1.99467e-2, $10*1.99467e-2}' >h2o.for
grep 'O\|H ' H2O.md |awk '(NR-1)%3==0 {printf "%.0f\n", NR/3+1} {print $2, $3, $4}' >h2o.pos

cat > h2o.in << EOF
 &control
 /
 &system
    ibrav = 14,
    celldm(1) = 12.0,
    celldm(2) = 1.0,
    celldm(3) = 1.0,
    celldm(4) = 0.0,
    celldm(5) = 0.0,
    celldm(6) = 0.0,
    nat  = 3,
    ntyp = 2,
 /
 &electrons
 /
 &ions
 /
 &cell
 /
ATOMIC_SPECIES
EOF
sed -n '/<Definition.of.Atomic.Species/,/Definition.of.Atomic.Species>/p' H2O.dat |sed '/</d; />/d' |awk '{print $1, "1.00d0", "H.blyp-vbc.UPF"}'>>h2o.in
echo "ATOMIC_POSITIONS (bohr)" >>h2o.in
sed -n '/<Atoms.SpeciesAndCoordinates/,/Atoms.SpeciesAndCoordinates>/p' H2O.dat |sed '/</d; />/d' |awk '{print $2, $3, $4, $5}'>>h2o.in