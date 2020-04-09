import annotate
import matplotlib as mpl
from matplotlib import cm
import pandas as pnd
import numpy as np


def example_ensemble():
    email = "janani.durairaj@gmail.com"
    # TODO: replace this with PDBe SIFT query
    pdb_ids = "6y2g, 5r7y, 5r7z, 5r80, 5r81, 5r82, 5r83, 5r84, 5re4, 5re5, 5re6, 5re7, 5re8, 5re9, 5rea, 5reb, 5rec, 5red, 5ree, 5ref, 5reg, 5reh, 5rei, 5rej, 5rek, 5rel, 5rem, 5ren, 5reo, 5rep, 5rer, 5res, 5ret, 5reu, 5rev, 5rew, 5rex, 5rey, 5rez, 5rf0, 5rf1, 5rf2, 5rf3, 5rf4, 5rf5, 5rf6, 5rf7, 5rf8, 5rf9, 5rfa, 5rfb, 5rfc, 5rfd, 5rfe, 5rff, 5rfg, 5rfh, 5rfi, 5rfj, 5rfk, 5rfl, 5rfm, 5rfn, 5rfo, 5rfp, 5rfq, 5rfr, 5rfs, 5rft, 5rfu, 5rfv, 5rfw, 5rfx, 5rfy, 5rfz, 5rg0, 6lu7, 6m03, 6w63, 6y2e, 6y2f, 6y84, 6yb7".split(", ")
    pdb_to_chain = {p: "A" for p in pdb_ids}
    
    for p in ["6y2e", "6y2f", "6y2g"]:
        pdb_to_chain[p] = "AAA" 
    pdb_chain_pairs = [(p, pdb_to_chain[p]) for p in pdb_ids]
    
    # Map PDB residue index to Reference residue (TODO: replace this with Genbank mapping query)
    uniprot_id = "P0DTD1"
    offset = 3264
    residue_mapper = dict(zip(range(306), range(offset, offset+306)))
    
    annotations = annotate.get_annotations_ensemble(uniprot_id, pdb_chain_pairs, residue_mapper)
    
    # Make post=True and change email to post to beta SWISS MODEL website
    annotators = annotate.make_swiss_model_annotators(annotations, post=True, email=email)
                          
    # This prints a lot, comment it out if posting
    for annotator in annotators:
         print(annotator)

            
def example_single():
    email = "janani.durairaj@gmail.com"
    pdb_id = "6y2g"
    chain = "AAA"
    
    # Map PDB residue index to Reference residue (TODO: replace this with Genbank mapping query)
    uniprot_id = "P0DTD1"
    offset = 3264
    residue_mapper = dict(zip(range(306), range(offset, offset+306)))
    annotations = annotate.get_annotations_single(uniprot_id, pdb_id, chain, residue_mapper, n_modes=6)
    
    # Make post=True and change email to post to beta SWISS MODEL website
    annotators = annotate.make_swiss_model_annotators(annotations, post=True, email=email)
                          
    # This prints a lot, comment it out if posting
    for annotator in annotators:
         print(annotator)
    

if __name__ == "__main__":
    example_single()
    example_ensemble()
    