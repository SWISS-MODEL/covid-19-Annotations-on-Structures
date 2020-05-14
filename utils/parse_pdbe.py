import ftplib
import gzip
import typing
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import intervaltree as it
import requests as rq

from utils.uniprot import seq_from_ac

MAPPING_FILE = "uniprot_segments_observed.tsv"


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
                    yield current_key, current_sequence
                    current_sequence = ""
                    current_key = line.split(">")[1].strip()
            else:
                current_sequence += line.strip()
        yield current_key, current_sequence


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


def get_sift_xml(pdb_id: str) -> ET.Element:
    """
    Gets XML file for a PDB ID from SIFTS via FTP

    Parameters
    ----------
    pdb_id

    Returns
    -------
    ElementTree parsed XML Element
    """
    ftp_address = "ftp.ebi.ac.uk"
    ftp = ftplib.FTP(ftp_address)
    ftp.login()
    filepath = f"/pub/databases/msd/sifts/split_xml/{pdb_id[1:3]}"
    filename = f"{pdb_id}.xml.gz"
    ftp.cwd(filepath)
    content = list()
    ftp.retrbinary(f"RETR {filename}", content.append)
    ftp.close()
    content = b''.join(content)
    return ET.fromstring(gzip.decompress(content).decode("utf-8"))


def get_pdb_to_uniprot_mapping(pdb_id: str):
    """
    Maps from PDB residue number to UniProt residue number for each chain
    Missing residues are ignored

    Parameters
    ----------
    pdb_id

    Returns
    -------
    dict of {chain: dict of {uniprot_resnum: pdb_resnum}}
    """
    sift_xml = get_sift_xml(pdb_id)
    entities = [x for x in sift_xml.iter() if "entity" in x.tag]
    chains = defaultdict(dict)
    for ent in entities:
        for residue in [x for x in ent.iter() if "residue" in x.tag and not x.tag.endswith("Detail")]:
            try:
                pdb_entry = [x for x in residue.iter() if "dbSource" in x.attrib and x.attrib["dbSource"] == "PDB"][0]
                pdb_resnum = pdb_entry.attrib["dbResNum"]
                pdb_chain_id = pdb_entry.attrib["dbChainId"]
                uniprot_resnum = [x for x in residue.iter() if "dbSource" in x.attrib and x.attrib["dbSource"] == "UniProt"][0].attrib["dbResNum"]
                if pdb_resnum != "null":
                    chains[pdb_chain_id][int(pdb_resnum)] = int(uniprot_resnum)
            except IndexError:
                continue
    return chains


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


def main():
    mapping = UniProtBasedMapping("P0DTC2")
    mapping.list_available_annotations()
    print(mapping.uniprot_sequence)
    print(len(mapping.uniprot_sequence))
    for x in mapping.search_pdbs_by_protein_name("Spike protein S1"):
        print(x)
    print(get_pdb_to_uniprot_mapping("6vsb"))


if __name__ == "__main__":
    main()
