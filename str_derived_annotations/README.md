# Extract and Visualize Structure-Derived Annotations from PDB Sets

## Library Dependencies

* `numpy`
* `prody`
* `matplotlib`
* `dataclasses`
* `dssp`
* `requests`

## Ensemble Based Output

For a set (ensemble) of related PDB IDs, e.g. the same protein crystallized with different ligands, the following annotations are produced is the directory set by the user:

* `RMSD_to_<PDB>_<CHAIN>.txt` : RMSD of each structure to the reference structure (the first PDB ID in list)
* `Average_RMSD.txt`: Average RMSD of each residue across all structures
* `PCA_fluctuations.txt`: Squared fluctuations of each residue according to a PCA across all structures

## Single PDB Based Output

For a single PDB ID, the following annotations are produced in the same directory. With the exception of solvent accessibility (see description below), elastic network models are utilized.

* `ENM_fluctuations.txt`: Residue fluctuations
* Perturbation response
  * `Perturbation_Effectiveness.txt`: Effectiveness of each residue in perturbing other residues
  * `Perturbation_Sensitivity.txt`: Sensitivity of each residue to perturbation
* `Mechanical_Stiffness.txt`: Mechanical stiffness
* `Hinge_sites_for_mode_X.txt`: Hinge sites connecting two stretches of structure that may move independently
* `Relative_Solvent_Accesibility.txt`: Relative solvent accessibility
* `Relative_Solvent_Accesibility_PDB_Selection.txt`: Description of PDB file and chain selection used in the calculation. The selection is particularly important for solvent accessibility (see below).

Solvent accessibility utilizes a single PDB and solvent exposed surface area predictions are as produced by [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html) and then normalized for each amino acid surface area, based on surface area as determined from the peptide G-X-G, as reported by [Chothia](https://www.sciencedirect.com/science/article/abs/pii/0022283676901911?via%3Dihub). It **IS** possible for relative surface area to exceed 1.0 using this method.

The solvent accessibility calculation is complicated in that it likely requires atoms beyond those identified by the reference mapper. For this reason, a separate variable (see `full_pdb_solvent_accessibility` in the `example_usage.py` script) is used for this calculation. When it is set to `True`, the entire PDB will be used. Otherwise (when it is `False`), only the chain(s) identified by the reference mapper are used for the calculation. Note that the output atoms are only those from the reference mapper, in all cases. The PDB ID and chain selection are noted in the output file `Relative_Solvent_Accesibility_PDB_Selection.txt`.

## Testing the Library

These annotations can be written out as structure annotation and visualized in the beta SWISS-MODEL annotation website (see `example_usage.py` and the annotations folder for examples).

