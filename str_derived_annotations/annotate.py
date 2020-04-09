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
from utils.sm_annotations import Annotation
import numpy as np
import prody as pd
from sys import argv
import matplotlib as mpl
from matplotlib import cm
from pathlib import Path
from dataclasses import dataclass
import typing


def get_structures(structure_chain_id_pairs, files_exist=False):
    pdb_to_chain = {p: c for p,c in structure_chain_id_pairs}
    if files_exist:
        return [pd.parseCIF(f"{x}.cif", chain=pdb_to_chain[x]) for x in pdb_to_chain.keys()]
    else:
        return [pd.parseCIF(x, chain=pdb_to_chain[x]) for x in pdb_to_chain.keys()]


def make_ensemble(structures):
    ensemble = pd.buildPDBEnsemble(structures, subset="calpha")
    ensemble.iterpose()
    return ensemble


def get_rmsds_to_reference(ensemble):
    return pd.calcRMSD(ensemble)


def get_pca_fluctuations(ensemble, limit=3):
    pca = pd.PCA()
    pca.buildCovariance(ensemble)
    pca.calcModes()
    return pd.calcSqFlucts(pca[:limit])


def get_enm_fluctuations(enm, n_modes=6):
    return pd.calcSqFlucts(enm[:n_modes])


def get_hinge_indices(enm, mode=3):
    return enm.getHinges(mode)


def get_rmsd_per_residue(ensemble):
    return ensemble.getMSFs()


def get_perturbations(enm, n_modes=6):
    _, effectiveness, sensitivity = pd.calcPerturbResponse(enm[:n_modes])
    return effectiveness, sensitivity


def get_stiffness(enm, calphas, n_modes=6):
    return np.mean(pd.calcMechStiff(enm[:n_modes], calphas), axis=0)


def numbers_to_colors(numbers, cmap="jet", log=False):
    if log:
        numbers = np.log1p(numbers)
    norm = mpl.colors.Normalize(vmin=np.min(numbers), vmax=np.max(numbers))
    colormap = cm.get_cmap(cmap)
    return [colormap(norm(n))[:3] for n in numbers]


@dataclass
class EnsembleAnnotation:
    pdb_id: str
    chain: str
    uniprot_id: str
    residue_mapper: dict
    rmsds_to_reference: typing.List[float]
    rmsds_per_residue: np.ndarray
    pca_fluctuations: np.ndarray
    ensemble: pd.PDBEnsemble
        
    def get_output_mapping(self):
        mapping = dict()
        msa = self.ensemble.getMSA().getArray()[0]
        reference_residues = [i for i, r in enumerate(msa) if r != "-"]
        mapped_residues = [self.residue_mapper[r] for r in reference_residues]
        mapping["Ensemble IDs"] = [(f"{mapped_residues[0]}-{mapped_residues[reference_residues[-1]]}", self.ensemble.getLabels(), (0, 0, 0))]
        mapping["Average RMSD"] = zip(mapped_residues,
                                      np.round(self.rmsds_per_residue[reference_residues], 2),
                                      numbers_to_colors(self.rmsds_per_residue[reference_residues]))
        mapping["PCA fluctuations"] = zip(mapped_residues,  
                                          np.round(self.pca_fluctuations[reference_residues], 2),
                                          numbers_to_colors(self.pca_fluctuations[reference_residues]))

        return mapping
    
@dataclass
class StructureAnnotation:
    pdb_id: str
    chain: str
    uniprot_id: str
    residue_mapper: dict
    enm_fluctuations: np.ndarray
    perturbation_effectiveness: np.ndarray
    perturbation_sensitivity: np.ndarray
    mechanical_stiffness: np.ndarray
    hinge_sites: list
    anm: pd.dynamics.anm.ANM
    gnm: pd.dynamics.gnm.GNM
        
    
    def get_output_mapping(self):
        mapping = dict()
        residues = [self.residue_mapper[i] for i in range(len(self.enm_fluctuations))]
        mapping["ENM fluctuations"] = zip(residues, np.round(self.enm_fluctuations, 2), numbers_to_colors(self.enm_fluctuations))
        mapping["Perturbation Effectiveness"] = zip(residues, np.round(self.perturbation_effectiveness, 2), numbers_to_colors(self.perturbation_effectiveness))
        mapping["Perturbation Sensitivity"] = zip(residues, np.round(self.perturbation_sensitivity, 2), numbers_to_colors(self.perturbation_sensitivity))
        mapping["Mechanical Stiffness"] = zip(residues, np.round(self.mechanical_stiffness, 2), numbers_to_colors(self.mechanical_stiffness))
        for i in range(len(self.hinge_sites)):
            mapping[f"Hinge sites for mode {i}"] = [(self.residue_mapper[r], f"mode {i}", (0, 0, 0) ) for r in self.hinge_sites[i]]
        return mapping

def get_annotations_ensemble(reference_uniprot_id, structure_chain_id_pairs, residue_mapper: list):
    structures = get_structures(structure_chain_id_pairs)
    ensemble = make_ensemble(structures)
    rmsds_to_reference = get_rmsds_to_reference(ensemble)
    rmsds_per_residue = get_rmsd_per_residue(ensemble)
    pca_fluctuations = get_pca_fluctuations(ensemble)
    return EnsembleAnnotation(structure_chain_id_pairs[0][0], structure_chain_id_pairs[0][1], reference_uniprot_id,
                              residue_mapper, rmsds_to_reference, rmsds_per_residue, pca_fluctuations, ensemble)


def get_annotations_single(uniprot_id, pdb_id, chain, residue_mapper, n_modes=6):
    structure = pd.parseCIF(pdb_id, chain=chain)
    gnm, calphas = pd.calcGNM(structure, n_modes=n_modes)
    anm, _ = pd.calcANM(structure, n_modes=n_modes)
    effectiveness, sensitivity = get_perturbations(anm, n_modes)
    hinge_sites = [get_hinge_indices(gnm, mode=n) for n in range(n_modes)]
    return StructureAnnotation(pdb_id, chain, uniprot_id, residue_mapper, get_enm_fluctuations(anm, n_modes), 
                               effectiveness, sensitivity, get_stiffness(anm, calphas, n_modes), hinge_sites, anm, gnm)


def make_swiss_model_annotators(annotation: typing.Union[StructureAnnotation, EnsembleAnnotation], 
                                post=False, email=None):
    output_mapping = annotation.get_output_mapping()
    annotators = []
    for title in output_mapping:
        annotator = Annotation()
        for residue, value, color in output_mapping[title]:
            annotator.add(annotation.uniprot_id, residue, tuple(color), str(value))
        annotators.append((title, annotator))
        if post:
            print(annotator.post(title=f"{title} PDB ID {annotation.pdb_id} Chain {annotation.chain}", email=email))
    return annotators