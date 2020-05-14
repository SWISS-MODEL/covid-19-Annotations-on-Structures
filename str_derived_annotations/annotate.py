# Plan:
# input:
# 1. ref protein id
# 2. list of related protein ids
# output per input protein and per residue:
# 1. signature fluctuations upto mode 6 (obtained from anm/gnm)
# 2. hinge sites
# 3. mean rmsd of each residue to the reference
# output per input protein
# 1. rmsd of each related protein to the reference protein
# output per reference protein and per residue
# 1. effectiveness of perturbation
# 2. sensitivity to perturbation
# 3. mechanical stiffness
# 4. pca fluctuations obtained from different confirmations
import os
import tempfile
import typing
from dataclasses import dataclass

import numpy as np
import prody as pd
from matplotlib import cm
from matplotlib import colors as mpl_colors

from utils.sm_annotations import Annotation


def get_structures(structure_chain_id_pairs: typing.List[typing.Tuple[str, str]]):
    """
    Gets ProDy AtomGroup objects for each (pdb_id, chain) pair
    """
    pdb_to_chain = {p: c for p, c in structure_chain_id_pairs}
    return [pd.parseCIF(x, chain=pdb_to_chain[x]) for x in pdb_to_chain.keys()]


def make_ensemble(structures: typing.List[pd.AtomGroup]):
    """
    Builds an Ensemble object that superposes all member structures onto each other
    """
    ensemble = pd.buildPDBEnsemble(structures, subset="calpha")
    ensemble.iterpose()
    return ensemble


def get_rmsds_to_reference(ensemble):
    """
    Gets RMSD of each structure to the reference
    """
    return pd.calcRMSD(ensemble)


def get_pca_fluctuations(ensemble, limit=3):
    """
    Get squared fluctuations of each residue according to a PCA on the ensemble
    Parameters
    ----------
    ensemble
        pd.PDBEnsemble object
    limit
        number of PCA modes to consider

    Returns
    -------
    array of squared fluctuations per aligned residue
    """
    pca = pd.PCA()
    pca.buildCovariance(ensemble)
    pca.calcModes()
    return pd.calcSqFlucts(pca[:limit])


def get_enm_fluctuations(enm, n_modes=6):
    """
    Get squared fluctuations of each residue according to an elastic network model
    Parameters
    ----------
    enm
        pd.dynamics.anm.ANM or pd.dynamics.gnm.GNM object
    n_modes
        number of ENM modes to consider

    Returns
    -------
    array of squared fluctuations per residue
    """
    return pd.calcSqFlucts(enm[:n_modes])


def get_hinge_indices(enm, mode):
    """
    Gets residue indices marked as hinge (based on ENM calculations for a specific mode)
    These are residues which act as a hinge between other (predicted) moving parts
    """
    return enm.getHinges(mode)


def get_rmsd_per_residue(ensemble):
    """
    Gets average RMSD of an aligned residue across proteins in an ensemble
    """
    return ensemble.getMSFs()


def get_perturbations(enm, n_modes=6):
    """
    Calculates perturbation response based on an elastic network model
    Parameters
    ----------
    enm
    n_modes
        number of modes to consider

    Returns
    -------
    Effectiveness - how effective each residue is in perturbing other residues
    Sensitivity - how sensitive a residue is to perturbation by other residues
    """
    _, effectiveness, sensitivity = pd.calcPerturbResponse(enm[:n_modes])
    return effectiveness, sensitivity


def get_stiffness(enm, calphas, n_modes=6):
    """
    Calculate mechanical stiffness based on an elastic network model
    Parameters
    ----------
    enm
    calphas
        alpha carbon selection from the structure
    n_modes
        number of modes to consider

    Returns
    -------
    a value of stiffness per residue
    """
    return np.mean(pd.calcMechStiff(enm[:n_modes], calphas), axis=0)


# Amino acid solvent accessibility measurements
# C. Chotia, The Nature of the Accessible and Buried Surfaces in Proteins, J. Mol. Biol., 105(1975)1-14.
# https://www.sciencedirect.com/science/article/abs/pii/0022283676901911?via%3Dihub
# http://prowl.rockefeller.edu/aainfo/volume.htm
AA_SA_VOL = dict(zip(['ALA', 'ARG', 'ASP', 'ASN', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                      'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
                     [ 115.,  225.,  150.,  160.,  135.,  190.,  180.,   75.,  195.,  175., 
                       170.,  200.,  185.,  210.,  145.,  115.,  140.,  255.,  230.,  155.]))

                               
def get_relative_solvent_accessibility(pdb_id, residue_mapper, chain, full_pdb_solvent_accessibility=True, aa_surface_area=AA_SA_VOL):
    """
    Run DSSP on a PDB file and return the resulting AtomGroup
    
    Parameters
    ----------
    pdb_id
        String containing PDB ID
    residue_mapper
        Dictionary of residue - unitprot mappings
    chain
        String containing the selected chain ID(s) from the residue mapper
    full_pdb_solvent_accessibility
        Boolean to use the full PDB for solvent accessibility calculations -- otherwise
        only the chain residues will be selected. Default is True.
    aa_surface_area
        Dictionary with amino acid abbreviations as keys and surface area 
        calculations as values

    Returns
    -------
    a numpy array containing relative solvent accessibility measurement for residues
    """
    
    if full_pdb_solvent_accessibility:
        dssp_chain = None
    else:
        dssp_chain = chain
        
    with tempfile.TemporaryDirectory() as tdir:
        pdb_file = os.path.join(tdir, '.'.join([pdb_id, 'pdb']))
        dssp_file = os.path.join(tdir, '.'.join([pdb_id, 'dssp']))
            
        # DSSP doesn't work with CIF-based atom groups, so must re-run here
        pd.pathPDBFolder(tdir)
        structure = pd.parsePDB(pdb_id, chain=dssp_chain)
        
        # Must write PDB file for DSSP with only chain selections
        # TODO how to silence output from the DSSP functions
        pd.writePDB(pdb_file, structure) 
        pd.execDSSP(pdb_file, outputdir=tdir) 
        pd.parseDSSP(dssp_file, structure)

    # Gather results
    # There should not be missing residues 
    mapped_residue_list = list(residue_mapper.keys())
    mapped_residue_list = ' '.join([str(x) for x in mapped_residue_list])

    selection_string = f"resnum {mapped_residue_list}"
    if dssp_chain is not None:
        selection_string += f" AND chain {chain}"
        
    iter_resi_list = sorted(set(structure.select(selection_string).getResnums()))
    rel_acc_list = list()

    for resi in iter_resi_list:
        dssp_resi = structure[(chain, resi)]
        surface_accessibilty = dssp_resi.getData('dssp_acc')[0]
        resn = dssp_resi.getResname()
        rel_surface_accessibilty = surface_accessibilty / aa_surface_area[resn]
        rel_acc_list.append(rel_surface_accessibilty)

    return np.array(rel_acc_list)


def numbers_to_colors(numbers, cmap="jet", log=False):
    """
    Converts a list of real-valued numbers to colors according to a colormap
    used for plotting on a structure in SWISS-MODEL
    Parameters
    ----------
    numbers
    cmap
        matplotlib colormap
    log
        if True, does natural log transformation on the numbers first
        (might be better if the distribution is too skewed)

    Returns
    -------

    """
    if log:
        numbers = np.log1p(numbers)
    norm = mpl_colors.Normalize(vmin=np.min(numbers), vmax=np.max(numbers))
    colormap = cm.get_cmap(cmap)
    return [colormap(norm(n))[:3] for n in numbers]


@dataclass
class EnsembleAnnotation:
    pdb_id: str
    chain: str
    protein: pd.AtomGroup

    uniprot_id: str
    residue_mapper: dict

    rmsds_to_reference: typing.List[float]
    rmsds_per_residue: np.ndarray
    pca_fluctuations: np.ndarray
    ensemble: pd.PDBEnsemble

    def get_output_mapping(self):
        mapping = dict()
        msa = self.ensemble.getMSA().getArray()[0]
        calphas = self.protein.select("calpha")
        reference_aln_indices = [i for i, r in enumerate(msa) if r != "-"]
        indices = []
        resnums = []
        for i, a in enumerate(calphas):
            if a.getResnum() in self.residue_mapper:
                indices.append(i)
                resnums.append(a.getResnum())
        uniprot_residues = [self.residue_mapper[r] for r in resnums]
        mapping["Ensemble IDs"] = [(f"{uniprot_residues[0]}-{uniprot_residues[-1]}",
                                    self.ensemble.getLabels(), (0, 0, 0))]
        mapping["Average RMSD"] = zip(uniprot_residues,
                                      np.round(self.rmsds_per_residue[reference_aln_indices], 2),
                                      numbers_to_colors(self.rmsds_per_residue[reference_aln_indices]))
        mapping["PCA fluctuations"] = zip(uniprot_residues,
                                          np.round(self.pca_fluctuations[reference_aln_indices], 2),
                                          numbers_to_colors(self.pca_fluctuations[reference_aln_indices]))

        return mapping

    def write_rmsds_to_reference(self, filename):
        with open(filename, "w") as f:
            f.write(f"PDB ID\tRMSD to reference PDB {self.pdb_id}_{self.chain}\n")
            for i, name in enumerate(self.ensemble.getLabels()[1:]):
                f.write(f"{name}\t{self.rmsds_to_reference[i + 1]}\n")


@dataclass
class StructureAnnotation:
    pdb_id: str
    chain: str
    protein: pd.AtomGroup
    calphas: pd.AtomGroup

    uniprot_id: str
    residue_mapper: dict

    enm_fluctuations: np.ndarray
    perturbation_effectiveness: np.ndarray
    perturbation_sensitivity: np.ndarray
    mechanical_stiffness: np.ndarray
    relative_solvent_accessibility: np.ndarray
    hinge_sites: list
    anm: pd.dynamics.anm.ANM
    gnm: pd.dynamics.gnm.GNM

    def get_output_mapping(self):
        mapping = dict()
        indices = []
        resnums = []
        for i, a in enumerate(self.calphas):
            if a.getResnum() in self.residue_mapper:
                indices.append(i)
                resnums.append(a.getResnum())
        uniprot_residues = [self.residue_mapper[i] for i in resnums]
        mapping["ENM fluctuations"] = zip(uniprot_residues, np.round(self.enm_fluctuations[indices], 2),
                                          numbers_to_colors(self.enm_fluctuations[indices]))
        mapping["Perturbation Effectiveness"] = zip(uniprot_residues, np.round(self.perturbation_effectiveness[indices], 2),
                                                    numbers_to_colors(self.perturbation_effectiveness[indices]))
        mapping["Perturbation Sensitivity"] = zip(uniprot_residues, np.round(self.perturbation_sensitivity[indices], 2),
                                                  numbers_to_colors(self.perturbation_sensitivity[indices]))
        mapping["Mechanical Stiffness"] = zip(uniprot_residues, np.round(self.mechanical_stiffness[indices], 2),
                                              numbers_to_colors(self.mechanical_stiffness[indices]))
        mapping["Relative Solvent Accessibility"] = zip(uniprot_residues, np.round(self.relative_solvent_accessibility[indices], 2),
                                                        numbers_to_colors(self.relative_solvent_accessibility[indices]))
        for i in range(len(self.hinge_sites)):
            mapping[f"Hinge sites for mode {i}"] = [(self.residue_mapper[self.calphas[r].getResnum()],
                                                     f"mode {i}",
                                                     (0, 0, 0)) for r in self.hinge_sites[i] if self.calphas[r].getResnum() in self.residue_mapper]
        return mapping


def get_annotations_ensemble(reference_uniprot_id, structure_chain_id_pairs, residue_mapper: dict):
    structures = get_structures(structure_chain_id_pairs)
    ensemble = make_ensemble(structures)
    rmsds_to_reference = get_rmsds_to_reference(ensemble)
    rmsds_per_residue = get_rmsd_per_residue(ensemble)
    pca_fluctuations = get_pca_fluctuations(ensemble)
    return EnsembleAnnotation(structure_chain_id_pairs[0][0], structure_chain_id_pairs[0][1], structures[0],
                              reference_uniprot_id, residue_mapper,
                              rmsds_to_reference, rmsds_per_residue, pca_fluctuations, ensemble)


def get_annotations_single(uniprot_id, pdb_id, residue_mapper: dict, chain=None, n_modes=6, full_pdb_solvent_accessibility=True):
    structure = pd.parseCIF(pdb_id, chain=chain)
    gnm, calphas = pd.calcGNM(structure, n_modes=n_modes)
    anm, _ = pd.calcANM(structure, n_modes=n_modes)
    effectiveness, sensitivity = get_perturbations(anm, n_modes)
    hinge_sites = [get_hinge_indices(gnm, mode=n) for n in range(n_modes)]
    return StructureAnnotation(pdb_id, chain, structure, calphas,
                               uniprot_id, residue_mapper,
                               get_enm_fluctuations(anm, n_modes), 
                               effectiveness, 
                               sensitivity, 
                               get_stiffness(anm, calphas, n_modes),
                               get_relative_solvent_accessibility(pdb_id, residue_mapper, chain=chain, 
                                                                  full_pdb_solvent_accessibility=full_pdb_solvent_accessibility),
                               hinge_sites, anm, gnm)


def make_swiss_model_annotators(annotation: typing.Union[StructureAnnotation, EnsembleAnnotation],
                                post=False, email=None):
    output_mapping = annotation.get_output_mapping()
    titles = []
    annotators = []
    urls = []
    for title in output_mapping:
        annotator = Annotation()
        for residue, value, color in output_mapping[title]:
            annotator.add(annotation.uniprot_id, residue, tuple(color), str(value))
        titles.append(title)
        annotators.append(annotator)
        if post:
            urls.append(annotator.post(title=f"{title} PDB ID {annotation.pdb_id} Chain {annotation.chain}", email=email))
    return annotators, titles, urls
