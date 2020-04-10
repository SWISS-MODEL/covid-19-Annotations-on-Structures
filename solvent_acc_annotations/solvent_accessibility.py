import os
from collections import defaultdict
from pathlib import Path
import numpy as np
import prody
import pandas as pd

# Amino acid solvent accessibility measurements
# http://prowl.rockefeller.edu/aainfo/volume.htm
AA_SA_VOL = pd.DataFrame({'resn':         ['ALA', 'ARG', 'ASP', 'ASN', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                                           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'],
                          'surface_area': [ 115.,  225.,  150.,  160.,  135.,  190.,  180.,   75.,  195.,  175., 
                                            170.,  200.,  185.,  210.,  145.,  115.,  140.,  255.,  230.,  155.],
                          'volume':       [ 88.6, 173.6, 111.1, 114.1, 108.5, 138.4, 143.8,  60.1, 153.2, 166.7, 
                                           166.7, 168.6, 162.9, 189.9, 112.7,  89.0, 116.1, 227.8, 193.6, 140.0]})


def run_dssp(pdb_id):
    """
    Run DSSP on a PDB file and return the resulting AtomGroup
    
    :param pdb_id: String containing PDB ID
    """
    # TODO: how to silence the C-level output from the DSSP functions completely
    
    file_ext_list = ['pdb', 'pdb.gz', 'dssp']
    gzip_file = pdb_id + '.pdb.gz'
    dssp_file = pdb_id + '.dssp'

    atom_group = prody.parsePDB(pdb_id, stderr=False)
    prody.execDSSP(gzip_file, stderr=False)
    prody.parseDSSP(dssp_file, atom_group)

    # File cleanup
    # TODO: is it possible to redirect file location in ProDyn? Would prefer to use tmpdir
    for ext in file_ext_list:
        filename = Path('.'.join([pdb_id, ext]))
        if filename.exists():
            filename.unlink()
            
    return atom_group


def gather_residue_dssp_data(atom_group):
    """
    Return DSSP data on a per-residue level and return as a dataframe
    
    :param atom_group: ProDyn AtomGroup instance
    """

    #attr_list = ['dssp_acc', 'dssp_kappa', 'dssp_alpha', 'dssp_phi', 'dssp_psi']
    attr_list = ['dssp_acc']
    dssp_dict = defaultdict(list)

    for res in atom_group.iterResidues():
        pdb_resi = res.getResindex()
        dssp_resi = int(np.mean(res.getData('dssp_resnum')))
        
        if dssp_resi != 0:
            dssp_dict['resi'].append(pdb_resi)
            dssp_dict['resn'].append(res.getResname())
            #dssp_dict['dssp_resi'].append(dssp_resi)
            
            for attr in attr_list:
                atom_data = np.array(res.getData(attr))
                assert np.allclose(atom_data.mean(), atom_data.min())
                dssp_dict[attr].append(atom_data.mean())

    dssp_df = pd.DataFrame(dssp_dict)
    return dssp_df


def calc_relative_solvent_sa(dssp_df, 
                             sa_df=AA_SA_VOL[['resn', 'surface_area']], 
                             resn_col='resn', resi_col='resi', sa_col='surface_area', rel_sa_col='dssp_acc_rel'):
    # TODO: is 
    """
    Calculate percentage surface area
    
    :param dssp_df: DataFrame with aggregated DSSP results -- output from `gather_residue_dssp_data`
    :param sa_df: DataFrame with per-residue surface accessible area measurements default is AA_SA_VOL
    :param resn_col: String that holds the column name for residue names
    :param resi_col: String thad holds the column name for residue numbers
    :param sa_col: String that holds the column name for surface area
    :param rel_sa_col: String that holds the column name 
    """
    rel_sa_df = pd.merge(dssp_df, sa_df, on=resn_col, how='left').sort_values(resi_col)
    rel_sa_df[rel_sa_col] = rel_sa_df['dssp_acc'] / rel_sa_df[sa_col]
    rel_sa_df.drop(sa_col, axis=1, inplace=True)
    return rel_sa_df


def write_results_to_csv(rel_sa_df, output_path):
    """
    Write calculation to CSV file
    
    :param rel_sa_df: DataFrame holding solvent accessibility calculations
    :param output_path: String holding the destination file name
    
    """
    
    column_dict = {'resi':'ResidueNumber', 'resn':'ResidueName', 
                 'dssp_acc':'SolventAcessibleArea', 'dssp_acc_rel':'RelativeAccessibleArea'}
    column_names = list(column_dict.values())
    
    rel_sa_df.rename(columns=column_dict)[column_names].to_csv(output_path, index=False)
    

def run_sa_calc(pdb_id, output_dir='.'):
    """
    Run the entire pipeline for a given PDB ID
    
    :param pdb_id: String containing PDB ID (no file extension)
    :param output_dir: String condaining output directory
    """
    atom_group = run_dssp(pdb_id)
    dssp_df = gather_residue_dssp_data(atom_group)
    rel_sa_df = calc_relative_solvent_sa(dssp_df)
    
    output_path = Path(os.path.join(output_dir, pdb_id + '_SA.csv'))
    write_results_to_csv(rel_sa_df, output_path)
    

if __name__ == '__main__':
    
    pdb_ids = "5r7y, 5r7z, 5r80, 5r81, 5r82, 5r83, 5r84, 5re4, 5re5, 5re6, 5re7, 5re8, 5re9, 5rea, 5reb, 5rec, 5red, 5ree, 5ref, 5reg, 5reh, 5rei, 5rej, 5rek, 5rel, 5rem, 5ren, 5reo, 5rep, 5rer, 5res, 5ret, 5reu, 5rev, 5rew, 5rex, 5rey, 5rez, 5rf0, 5rf1, 5rf2, 5rf3, 5rf4, 5rf5, 5rf6, 5rf7, 5rf8, 5rf9, 5rfa, 5rfb, 5rfc, 5rfd, 5rfe, 5rff, 5rfg, 5rfh, 5rfi, 5rfj, 5rfk, 5rfl, 5rfm, 5rfn, 5rfo, 5rfp, 5rfq, 5rfr, 5rfs, 5rft, 5rfu, 5rfv, 5rfw, 5rfx, 5rfy, 5rfz, 5rg0, 6lu7, 6m03, 6w63, 6y2e, 6y2f, 6y2g, 6y84, 6yb7"
    pdb_ids = pdb_ids.split(", ")
    
    output_dir = './output'
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    for pdb_id in pdb_ids:
        run_sa_calc(pdb_id, output_dir)
