from pathlib import Path
import numpy as np
from str_derived_annotations import annotate, parse_pdbe


def create_csv_file(annotations, csv_dir):
    """
    Create CSV file for annotations -- currently this is just implemented single PDB annotations
    
    Parameters
    ----------
    annotations
        Annotations object
    csv_dir
        Path containing output directory for CSV files
    """ 
    
    keys = ['enm_fluctuations', 'perturbation_effectiveness', 'perturbation_sensitivity', 
            'mechanical_stiffness', 'relative_solvent_accessibility']
        
    def all_same_len(items):
        return all([len(x)==len(items[0]) for x in items])

    data_dict = vars(annotations)
    data = {'residue_numbers': list(data_dict['residue_mapper'].keys()),
            'uniprot_annotation': list(data_dict['residue_mapper'].values())}
    for key in keys:
        data[key] = data_dict[key]

    assert all_same_len(list(data.values()))
    n_rows = len(data[keys[0]])

    col_names = list(data.keys())
    nice_title = "_".join([data_dict['uniprot_id'], data_dict['pdb_id']])
    if 'chain' in data_dict:
        nice_title += "_".join(['', 'chain', data_dict['chain']])
        
    with open(csv_dir / f"{nice_title}.csv", "w") as f:
        f.write(','.join(col_names) + "\n")
        for row in range(n_rows):
            data_row = [np.round(data[col][row],2) for col in col_names]
            f.write(','.join([str(x) for x in data_row]) + '\n')
            

def create_annotation_file(annotators, titles, urls, annotation_dir):
    """
    Create annotation files that can be uploaded to SWISS-MODEL
    
    Parameters
    ----------
    annotators
        Object created by the annotation class which contains the annotaiton values
    titles
        List of string titles for the annotations, used for filenames
    urls
        List of URL strings containing information (PDB ID, notes)
    annotation_dir
        Directory that SWISS-MODEL files will be written to
    """
    
    for i, (title, annotator) in enumerate(zip(titles, annotators)):
        nice_title = '_'.join(title.split())
        with open(annotation_dir / f"{nice_title}.txt", "w") as f:
            if len(urls):
                # This writes the URL to the first line of the file.
                # TODO: check if # can be used as a comment in the SWISS-MODEL annotation file
                f.write(f"# URL: {urls[i]}\n")
            f.write(str(annotator))
            
                
def example_ensemble(uniprot_id, pdb_id, pdb_search_str, output_path, post=False, email=None):
    """
    Create ensemble, calculate annotations, and upload to SWISS MODEL
    
    Parameters
    ----------
    uniprot_id
        String containing uniprot id
    pdb_id
        String containing reference PDB ID
    pdb_search_str
        String that will be used to search for other structures to include in the ensemble
    output_path
        String containing the output directory
    post
        Boolean that controls whether results are sent to SWISS MODEL
        Default is False, setting to True requires email
    email
        String containing email or None, required to post results
        Default is None
    """

    if post:
        assert email is not None, 'Email must be set to post results'
        
    reference_mapper = parse_pdbe.UniProtBasedMapping(uniprot_id)
    pdb_info_list = reference_mapper.search_pdbs_by_protein_name(pdb_search_str)
    pdb_chain_pairs = [(p["pdb_id"], p["chain_id"]) for p in pdb_info_list]
    residue_mapping, _ = reference_mapper.map_to_pdb(pdb_info_list[0])
    annotations = annotate.get_annotations_ensemble(uniprot_id, pdb_chain_pairs, residue_mapping)

    # Create output directory
    annotation_dir = Path(output_path)
    if not annotation_dir.exists():
        annotation_dir.mkdir(parents=True)

    # Make post=True and change email to post to beta SWISS MODEL website
    annotators, titles, urls = annotate.make_swiss_model_annotators(annotations, post=post, email=email)
    create_annotation_file(annotators, titles, urls, annotation_dir)


def example_single(pdb_id, uniprot_id, output_path, post=False, email=None):
    """
    Calculate annotations for single PDB and upload to SWISS MODEL
    
    Parameters
    ----------
    pdb_id
        String containing PDB ID
    uniprot_id
        String containing uniprot id
    output_path
        String containing the output directory
    post
        Boolean that controls whether results are sent to SWISS MODEL
        Default is False, setting to True requires email
    email
        String containing email or None, required to post results
        Default is None
    """
    # TODO fix  / test chain=None
        
    if post:
        assert email is not None, 'Email must be suplied to post results'
    
    reference_mapper = parse_pdbe.UniProtBasedMapping(uniprot_id)
    pdb_info = reference_mapper.search_pdb_by_id(pdb_id)
    residue_mapping, _ = reference_mapper.map_to_pdb(pdb_info)
    annotations = annotate.get_annotations_single(uniprot_id, pdb_id, residue_mapping, pdb_info["chain_id"], n_modes=6)

    # Create output directory
    annotation_dir = Path(output_path)
    if not annotation_dir.exists():
        annotation_dir.mkdir(parents=True)
        
    # Write output as CSV file
    create_csv_file(annotations, annotation_dir)
    
    # Make post=True and change email to post to beta SWISS MODEL website
    annotators, titles, urls = annotate.make_swiss_model_annotators(annotations, post=post, email=email)
    create_annotation_file(annotators, titles, urls, annotation_dir)
    # print(list(zip(titles, urls)))


if __name__ == "__main__":
    pdb_id = "6m71"
    uniprot_id = "P0DTD1"
    pdb_search_str = "3C-like proteinase"
    output_path = "str_derived_annotations/annotations"
    post = False
    email = None
    
    example_single(pdb_id, uniprot_id, output_path, post, email)
    example_ensemble(uniprot_id, pdb_id, pdb_search_str, output_path, post, email)
