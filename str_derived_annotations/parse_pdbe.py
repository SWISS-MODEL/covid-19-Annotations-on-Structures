import subprocess
import typing
from pathlib import Path

import intervaltree as it
import prody as pd
import requests as rq

from utils.uniprot import seq_from_ac


def get_sequences_from_fasta_yield(fasta_file: typing.Union[str, Path]) -> tuple:
    """
    Returns (accession, sequence) iterator
    Parameters
    ----------
    fasta_file
    Returns
    -------
    (accession, sequence)
    """
    with open(fasta_file) as f:
        current_sequence = ""
        current_key = None
        for line in f:
            if not len(line.strip()):
                continue
            if "==" in line:
                continue
            if ">" in line:
                if current_key is None:
                    current_key = line.split(">")[1].strip()
                else:
                    yield (current_key, current_sequence)
                    current_sequence = ""
                    current_key = line.split(">")[1].strip()
            else:
                current_sequence += line.strip()
        yield (current_key, current_sequence)


def get_sequences_from_fasta(fasta_file: typing.Union[str, Path]) -> dict:
    """
    Returns dict of accession to sequence from fasta file
    Parameters
    ----------
    fasta_file
    Returns
    -------
    {accession:sequence}
    """
    return {key: sequence for (key, sequence) in get_sequences_from_fasta_yield(fasta_file)}


def run_muscle(in_file, out_file):
    command = ['muscle', '-in', in_file, '-out', out_file]
    subprocess.check_call(command)


class UniProtBasedMapping:
    def __init__(self, uniprot_id: str):
        self.uniprot_id = uniprot_id
        self.PDBe_api_request_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
        self.uniprot_api_request_url = f"https://www.ebi.ac.uk/uniprot/api/covid-19/uniprotkb/accession/{uniprot_id}.gff"
        self.uniprot_sequence = seq_from_ac(uniprot_id)
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
            raise KeyError(f"No such protein name found. here are the available ones: \
            {', '.join(list(self.protein_annotation_intervals.keys()))}")
        return [self.data[x.data] for x in self.tree.overlap(start, end)]

    def search_pdb_by_id(self, pdb_id):
        try:
            return self.data[pdb_id]
        except KeyError:
            raise KeyError(f"PDB ID {pdb_id} not found. "
                           f"Are you sure this is a COVID-19 protein present in {self.uniprot_id}? "
                           f"Available IDs are {self.data.keys()}")

    def map_to_pdb(self, info_dict, aln_dir="./tmp"):
        pdb_id, chain_id = info_dict["pdb_id"], info_dict["chain_id"]
        pdb = pd.parseCIF(pdb_id, chain=chain_id).select(f"resnum {info_dict['start']}:{info_dict['end']}")
        pdb_alpha = pdb.select("calpha")
        pdb_sequence = pdb_alpha.getSequence()
        ref_sequence = self.uniprot_sequence[info_dict["unp_start"] - 1: info_dict["unp_end"]]
        aln_dir = Path(aln_dir)
        if not aln_dir.exists():
            aln_dir.mkdir()
        sequence_file = aln_dir / f"{pdb_id}_{chain_id}.fasta"
        with open(sequence_file, "w") as f:
            f.write(f">{pdb_id}\n{pdb_sequence}\n")
            f.write(f">{self.uniprot_id}\n{ref_sequence}\n")
        aln_file = aln_dir / f"{pdb_id}_{chain_id}_aln.fasta"
        run_muscle(sequence_file, aln_file)
        aln_sequences = get_sequences_from_fasta(aln_file)
        index_1 = info_dict["start"]
        index_2 = info_dict["unp_start"]
        residue_mapping = {}
        for i in range(len(aln_sequences[pdb_id])):
            if aln_sequences[pdb_id][i] != "-":
                residue_mapping[index_1] = index_2
                index_1 += 1
            if aln_sequences[self.uniprot_id] != "-":
                index_2 += 1
        return residue_mapping


if __name__ == "__main__":
    mapping = UniProtBasedMapping("P0DTC2")
    mapping.list_available_annotations()
    print(mapping.uniprot_sequence)
    print(len(mapping.uniprot_sequence))
    print(mapping.search_pdbs_by_protein_name("Spike protein S1"))
    # mapping.protein_annotation_intervals["3C-like proteinase"]
    # print(mapping.data["6vxx"])
