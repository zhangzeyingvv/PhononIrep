# PhononIrep

A phonon irreducible representations calculator([arXiv2201.11350](https://arxiv.org/pdf/2201.11350.pdf)).

# Installation
One need install both Python and Mathematica.
## Configure Python
First you should install python3 and phonopy.
e.g. install anaconda or miniconda, and run 
```
conda install -c conda-forge phonopy
```
in the anaconda(miniconda) prompt. 

To interface with Mathematica one need also install pyzmq by
```
pip install pyzmq 
```
## Configure Mathematica
According to the test, PhonoIrep cannot work in Matnematica 12.0, Therefore I recommended to run PhonoIrep in Mathematica 12.2 or 12.3. I will try to fix compatibility issues in the furture. 

Copy the SpaceGroupIrep direcory and PhononIrep.wl into any of the Mathematica `$Path` direcory.

That all for installation.



# Run PhononIrep


You can open a Mathematica notebook, save it in the directory where you do your phonopy calculations, and run
```
path = "D:\\Anaconda3"
SetEnvironment["PATH" -> Environment["PATH"] <> ";" <> path <> "\\Library\\bin"];
RegisterExternalEvaluator["Python", path <> "\\python.exe"];
Needs["PhononIrep`"]
```
to import PhononIrep. Here, you need to modify the `path` to the Python installation direcory. For Linux and Mac OS `python.exe` should be modfied by `python`, and run 
```
calcPhononIrep["supercell" -> {2, 2, 1},
 "unitcell" -> NotebookDirectory[] <> "POSCAR-unitcell",
 "force" ->  NotebookDirectory[] <> "FORCE_CONSTANTS",
 "kset" -> {{0,0,0}},
 "showRep" -> True
 ]
 ```
 to calculate the irreducible representations of Gamma point for input structure.
 
# Release Notes

v1.00b 2022/06/17

Add an example for diamond.

