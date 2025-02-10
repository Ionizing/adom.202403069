#!/bin/bash -e

for i in {001..018}; do
    echo $i;
    cp POSCAR-${i}  disp-${i}/POSCAR
    cp KPOINTS POTCAR INCAR sub_vasp_6.3 disp-${i}
done
