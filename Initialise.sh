#!/bin/bash

# Exit immediately if a command fails
set -e

nodes=$(lua -e "dofile('parameters.lua'); print(nodes)")
echo "Running on $nodes nodes"
echo "----------------------------------"
echo "Starting Initialisation"
echo ""
module load elmerfem


echo ""
echo "Building Mesh"

python3 Initialisation/createGeoFile.py
gmsh footprint.geo -1 -2 -o footprint.msh -v 0

ElmerGrid 14 2 footprint.msh -partition $nodes 1 1 -parttol 1 -autoclean  > /dev/null 2>&1
ElmerGrid 14 5 footprint.msh -partition $nodes 1 1 -parttol 1 -autoclean  > /dev/null 2>&1

echo "Done"
echo""

SCRIPT_FOLDER="/import/ontap-m-glaciology/Lloyd/Flowline/FlowlineScaled/Initialisation"
cd "$SCRIPT_FOLDER" || { echo "Folder $SCRIPT_FOLDER not found"; exit 1; }


echo "Scaling Initial Conditions"
python3 InitialiseGeometry.py
echo "Done"
echo ""
cd ../

#echo "Solving for Steady-State Strain Heating"
#nohup mpirun -np $nodes ElmerSolver InitialiseVel.sif > Initialisation/ScriptOutputs/InitialVel.out 
#echo "Done"
#echo ""

#echo "Solving for Diffusion Only Temperature Field"
#nohup mpirun -np $nodes ElmerSolver InitialiseT1.sif > Initialisation/ScriptOutputs/InitialT1.out 
#echo "Done"
#echo ""

cd SRC
elmerf90 sheetsolverhw.f90 -o SheetSolverhw.so # > /dev/null 2>&1
cd ../

echo "Solving for Full Temperature Field"
nohup mpirun -np $nodes ElmerSolver InitialiseT2.sif > Initialisation/ScriptOutputs/InitialT2.out 
echo "Done"
echo ""
echo "Initial Conditions Calculated"
