import numpy as np
import intervaltree as it
import requests as rq


class UniprotBasedMapping:
    def __init__(self, uniprot_id: str):
        self.PDBe_api_request_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
        self.uniprot_api_request_url = f"https://www.ebi.ac.uk/uniprot/api/covid-19/uniprotkb/accession/{uniprot_id}.gff"
        self.protein_annotation_intervals = dict()
        self._get_intervals_from_uniprot()
        self.data = {x["pdb_id"]: x for x in rq.get(self.PDBe_api_request_url).json()[uniprot_id]}
        self.tree = it.IntervalTree()
        self._build_tree()

    def list_available_annotations(self):
        print('\n'.join(list(self.protein_annotation_intervals.keys())))

    def _build_tree(self):
        for pdb_id, data in self.data.items():
            self.tree[data["unp_start"]: data["unp_end"]] = pdb_id

    def _get_intervals_from_uniprot(self):
        gff_lines = [x for x in rq.get(self.uniprot_api_request_url).text.split("\n") if not x.startswith("#") and len(x)]
        for line in gff_lines:
            line = line.split("\t")
            _range = (int(line[3]), int(line[4]))
            info = line[-2].split(";")
            note = [x for x in info if x.startswith("Note")]
            if len(note):
                protein_id = note[0].split("=")[-1]
                self.protein_annotation_intervals[protein_id] = _range
            else:
                continue

    def search_pdbs_by_protein_name(self, annotation_name):
        try:
            start, end = self.protein_annotation_intervals[annotation_name]
        except KeyError:
            raise KeyError(f"no such protein name found. here are the available ones: \
            {', '.join(list(self.protein_annotation_intervals.keys()))}")
        return [x.data for x in self.tree.overlap(start, end)]


if __name__ == "__main__":
    mapping = UniprotBasedMapping("P0DTC2")
    mapping.list_available_annotations()
    print(mapping.search_pdbs_by_protein_name("Spike protein S1"))
    # mapping.protein_annotation_intervals["3C-like proteinase"]