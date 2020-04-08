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

import numpy as np
import prody as pd
from sys import argv
import matplotlib as mpl
from pathlib import Path


def get_structures(structure_chain_id_pairs, files_exist=False):
    pdb_to_chain = {p: c for p,c in structure_chain_id_pairs}
    if files_exist:
        return [pd.parseCIF(f"{x}.cif", chain=pdb_to_chain[p]) for x in pdb_to_chain.keys()]
    else:
        return [pd.parseCIF(x, chain=pdb_to_chain[p]) for x in pdb_to_chain.keys()]


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

def get_enms(ensemble, method="anm"):
    return pd.calcEnsembleENMs(ensemble, model=method).getModeSets()

def get_enm_fluctuations(enms, n_modes=6):
    return [pd.calcSqFlucts(enms[i][:n_modes]) for i in range(len(enms))]


def get_hinge_indices(enms, mode=3):
    return [enms[i].getHinges(mode) for i in range(ensemble.numConfs())]


def get_rmsd_per_residue(ensemble):
    indices = ensemble._indices
    if indices is None:
        coords = ensemble._coords
        confs = ensemble._confs
        weights = ensemble._weights > 0
    else:
        coords = ensemble._coords[indices]
        confs = ensemble._confs[:, indices]
        weights = ensemble._weights[:, indices] > 0
    weightsum = weights.sum(0)
    mean = np.zeros(coords.shape)
    for i, conf in enumerate(confs):
        mean += conf * weights[i]
    mean /= weightsum

    ssqf = np.zeros(mean.shape)
    stds = []
    for i, conf in enumerate(confs):
        ssqf += ((conf - mean) * weights[i]) ** 2
        stds.append(((conf - mean) * weights[i]) ** 2)
    return ssqf.sum(1) / weightsum.flatten(), np.std(np.array(stds), axis=0).sum(1)


def get_perturbations(enms, n_modes=6):
    effectiveness = []
    sensitivity = []
    for i in range(len(enms)):
        _, e, s = pd.calcPerturbResponse(enms[i][:n_modes])
        effectiveness.append(e)
        sensitivity.append(s)
    return effectiveness, sensitivity


def get_stiffness(enms, structures, n_modes=6):
    return [np.mean(pd.calcMechStiff(enms[i][:n_modes], structures[i].select("calpha")), axis=0) for i in range(len(enms))]


def numbers_to_colors(numbers, cmap="jet", log=False):
    if log:
        numbers = np.log1p(numbers)
    norm = mpl.colors.Normalize(vmin=np.min(numbers), vmax=np.max(numbers))
    colormap = mpl.cm.get_cmap(cmap)
    return [colormap(norm(n))[:3] for n in numbers]


def main(structure_chain_id_pairs, output_dir):
    structures = get_structures(structure_chain_id_pairs)
    ensemble = make_ensemble(structures)
    rmsds_to_reference = get_rmsds_to_reference(ensemble)
    rmsds_per_residue, _ = get_rmsd_per_residue(ensemble)
    pca_fluctuations = get_pca_fluctuations(ensemble)

    enms = get_enms(ensemble)
    enm_fluctuations = get_enm_fluctuations(enms)
    effectiveness, sensitivity = get_perturbations(enms)
    stiffness = get_stiffness(enms, structures)

    output_dir = Path(output_dir)

    with open(output_dir / f"{'_'.join(structure_chain_id_pairs[0])}_reference.txt", "w") as f:
        f.write("Aligned Index\tResidue Index\tRef AA\tRMSDs per residue\tPCA fluctuations\n")
        msa = ensemble.getMSA().getArray()[0]
        j = 0
        for i in range(len(rmsds_per_residue)):
            f.write("\t".join(str(x) for x in [i] + [j, msa.astype("U")[i], rmsds_per_residue[i], pca_fluctuations[j]]) + "\n")
            if msa[i] != "-":
                j += 1
    for n in range(len(structure_chain_id_pairs)):
        with open(output_dir / f"{'_'.join(structure_chain_id_pairs[n])}.txt", "w") as f:
            f.write("Index\tANM fluctuations\tEffectiveness\tSensitivity\tStiffness\n")
            for i in range(len(enm_fluctuations[n])):
                f.write("\t".join(str(x) for x in [i] + [enm_fluctuations[n][i], effectiveness[n][i], sensitivity[n][i], stiffness[n][i]]) + "\n")


if __name__ == "__main__":
    # first ID is the reference
    pdb_ids = "6y2g, 5r7y, 5r7z, 5r80, 5r81, 5r82, 5r83, 5r84, 5re4, 5re5, 5re6, 5re7, 5re8, 5re9, 5rea, 5reb, 5rec, 5red, 5ree, 5ref, 5reg, 5reh, 5rei, 5rej, 5rek, 5rel, 5rem, 5ren, 5reo, 5rep, 5rer, 5res, 5ret, 5reu, 5rev, 5rew, 5rex, 5rey, 5rez, 5rf0, 5rf1, 5rf2, 5rf3, 5rf4, 5rf5, 5rf6, 5rf7, 5rf8, 5rf9, 5rfa, 5rfb, 5rfc, 5rfd, 5rfe, 5rff, 5rfg, 5rfh, 5rfi, 5rfj, 5rfk, 5rfl, 5rfm, 5rfn, 5rfo, 5rfp, 5rfq, 5rfr, 5rfs, 5rft, 5rfu, 5rfv, 5rfw, 5rfx, 5rfy, 5rfz, 5rg0, 6lu7, 6m03, 6w63, 6y2e, 6y2f, 6y84, 6yb7".split(", ")
    pdb_to_chain = {p: "A" for p in pdb_ids}
    for p in ["6y2e", "6y2f", "6y2g"]:
        pdb_to_chain[p] = "AAA"
    pdb_chain_pairs = [(p, pdb_to_chain[p]) for p in pdb_ids]
    main(pdb_chain_pairs, "annotations")
