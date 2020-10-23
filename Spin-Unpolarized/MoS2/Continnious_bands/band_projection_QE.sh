#!/bin/bash

#--------------------------USER DEFINED---------------
MPI="aprun -n"
CORES=40
#-------------------END OF USER DEFINITION------------

MPIRUN="${MPI} ${CORES}"

$MPIRUN  pw.x < scf.in | tee scf.out

$MPIRUN  pw.x < bands.in | tee bands.out

$MPIRUN  bands.x < bandsx.in | tee bandsx.out

$MPIRUN  projwfc.x < projwfc.in | tee projwfc.out
