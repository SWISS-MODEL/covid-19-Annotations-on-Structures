import re
import urllib.request

"""
Collection of handy functions related to uniprot. Potential reimplementations
of code that would be available in various packages with the goal of keeping
dependencies at a minimum.
"""


def valid_uniprot_ac_pattern(uniprot_ac):
    """
    Checks whether Uniprot AC is formally correct according to
    https://www.uniprot.org/help/accession_numbers 
    This is no check whether it actually exists.

    :param uniprot_ac:  Accession code to be checked

    """
    ac_pat = "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
    if re.match(ac_pat, uniprot_ac):
        return True
    else:
        return False


def seq_from_ac(uniprot_ac):
    """
    Fetches raw sequence string for given uniprot accession code

    :param uniprot_ac:  Accession code for which you want the sequence
    """
    if not valid_uniprot_ac_pattern(uniprot_ac):
        raise RuntimeError("Uniprot AC does not look valid")

    data = None
    try:
        # that's the default uniprot access
        url = "https://www.uniprot.org/uniprot/%s.fasta" % uniprot_ac
        with urllib.request.urlopen(url) as response:
            data = response.readlines()

    except:
        # this is only temporary, as SARS-CoV2 is not yet in uniprot
        url = (
            "https://www.ebi.ac.uk/uniprot/api/covid-19/uniprotkb/accession/%s.fasta"
            % (uniprot_ac)
        )
        with urllib.request.urlopen(url) as response:
            data = response.readlines()

    return "".join(line.decode().strip() for line in data[1:])
