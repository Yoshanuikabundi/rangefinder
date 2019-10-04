# rangefinder
Rangefinder scores residues of a provided solute binding protein (SBP) to help 
design semi-synthetic FRET sensors. 

Rangefinder depends on Python 3 (https://www.python.org/) and Biopython 
(http://biopython.org/wiki/Documentation), which should be installed before 
Rangefinder is used.

Rangefinder requires reliable structural models of the target SBP in both its 
open and closed conformation.

Rangefinder places a hypothetical fluorophore bead a set distance (default 
20 Angstroms) from the crystal structure's first CA atom, such that the bead,
the first CA atom, and the structure's centre of geometry form a line with the
CA atom in the centre. It then fits the two structures by the residues passed 
to the --align argument, and predicts dynamic ranges for FRET constructs with 
an ECFP fluorophore at the bead and an Alexa Fluor 532 dye at the CA atom of 
each residue.

First, acquire models of your SBP in PDB format in both conformations, and 
ensure their residue numbering is consistent. Then determine the residues to 
align the structures by. For sensors constructed by fusing a fluorescent 
protein to the N-terminus of a SBP, this should be the N-terminal lobe of the 
SBP.

Call Rangefinder with the following:

```python3 rangefinder.py -a <open structure> -b <closed structure> -f <residues to fit>```

Additional help and arguments are available with:

```python3 rangefinder.py --help```

As an example, maltose binding protein structures and output are provided. 
`example/1anf.pdb` and `example/1jw4.pdb` are taken from the PDB; 
`example/sensor_predict.csv` was produced by the command:

```python3 ../rangefinder.py -a 1jw4.pdb -b 1anf.pdb -f 1-110,259-285```
