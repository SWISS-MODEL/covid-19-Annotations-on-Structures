from pathlib import Path

import annotate
import parse_pdbe


def example_ensemble():
    uniprot_id = "P0DTD1"
    reference_mapper = parse_pdbe.UniProtBasedMapping(uniprot_id)
    pdb_info_list = reference_mapper.search_pdbs_by_protein_name("3C-like proteinase")
    pdb_chain_pairs = [(p["pdb_id"], p["chain_id"]) for p in pdb_info_list]
    residue_mapping = reference_mapper.map_to_pdb(pdb_info_list[0])
    annotations = annotate.get_annotations_ensemble(uniprot_id, pdb_chain_pairs, residue_mapping)

    # Make post=True and change email to post to beta SWISS MODEL website
    annotators, titles, urls = annotate.make_swiss_model_annotators(annotations, post=False, email=None)
    annotation_dir = Path("annotations")
    if not annotation_dir.exists:
        annotation_dir.mkdir()
    for i, (title, annotator) in enumerate(zip(titles, annotators)):
        nice_title = '_'.join(title.split())
        with open(annotation_dir / f"{nice_title}.txt", "w") as f:
            if len(urls):
                # This writes the URL to the first line of the file.
                # TODO: check if # can be used as a comment in the SWISS-MODEL annotation file
                f.write(f"# URL: {urls[i]}\n")
            f.write(str(annotator))


def example_single():
    pdb_id = "6m71"
    uniprot_id = "P0DTD1"
    reference_mapper = parse_pdbe.UniProtBasedMapping(uniprot_id)
    pdb_info = reference_mapper.search_pdb_by_id(pdb_id)
    residue_mapping = reference_mapper.map_to_pdb(pdb_info)
    annotations = annotate.get_annotations_single(uniprot_id, pdb_id, pdb_info["chain_id"], residue_mapping, n_modes=6)

    # Make post=True and change email to post to beta SWISS MODEL website
    annotators, titles, urls = annotate.make_swiss_model_annotators(annotations, post=False, email=None)
    annotation_dir = Path("annotations")
    if not annotation_dir.exists:
        annotation_dir.mkdir()
    for i, (title, annotator) in enumerate(zip(titles, annotators)):
        nice_title = '_'.join(title.split())
        with open(annotation_dir / f"{nice_title}.txt", "w") as f:
            if len(urls):
                # This writes the URL to the first line of the file.
                # TODO: check if # can be used as a comment in the SWISS-MODEL annotation file
                f.write(f"# URL: {urls[i]}\n")
            f.write(str(annotator))
    # print(list(zip(titles, urls)))


if __name__ == "__main__":
    example_single()
    example_ensemble()
