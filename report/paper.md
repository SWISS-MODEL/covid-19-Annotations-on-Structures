---
title: 'Mapping sequence data onto protein structures of SARS-CoV-2'
title_short: 'Mapping sequence data onto protein structures of SARS-CoV-2'
tags:
  - Covid-19
  - Protein structures
  - Sequence annotations
authors:
  - name: Gerardo Tauriello
    orcid: 0000-0002-5921-7007
    affiliation: 1
  - name: Xavier Robin
    affiliation: 1
  - name: Stefan Bienert
    affiliation: 1
  - name: Andrew Waterhouse
    affiliation: 1
  - name: Gabriel Studer
    affiliation: 1  
  - name: Laura Blum
    affiliation: 2
  - name: Michelle Gill
    orcid: 0000-0002-2861-355X
    affiliation: 3
    
affiliations:
 - name: SIB Swiss Institute of Bioinformatics & Biozentrum, University of Basel, Klingelbergstrasse 50–70, CH-4056 Basel, Switzerland
   index: 1
 - name: Middlebury College, Middlebury, VT, 05753, United States of America
   index: 2
 - name: NVIDIA, 2788 San Tomas Expressway, Santa Clara, CA, 95051, United States of America
   index: 3
date: 11 April 2020
bibliography: paper.bib
event: Covid2020
authors_short: Gerardo Tauriello \emph{et al.}
---

# Introduction

TODO:
- Description of problem
- Quick intro to annotation system based on SWISS-MODEL Repository [@Bienert-2017; @Waterhouse-2018]
- Introduce work done during BioHackathon

# Hackathon results

TODO:
- Subsection for each annotation effort

## Collection of annotation ideas

## Example code for annotations

To showcase the effort and provide an example of processing structures with
OpenStructure [@Biasini-2013], we extracted nucleotide binding residues for the
RNA-Polymerase as well as differences between SARS-CoV and SARS-CoV-2 in that
protein. This allowed to display both variations and the binding site on a
structure of the RNA-Polymerase (see figure \ref{RdRp_combined}). The variations
are extracted by aligning the UniProt-sequences P0DTD1 and P0C6X7 with a
pairwise Needleman-Wunsch algorithm and extracting an annotation for each
non-conserved amino acid. Insertions and deletions were ignored. An EM
structure exists for the RNA polymerase (PDB ID 6m71) but it does not contain
any ligands. We use structures from a study of viral RNA-Polymerases in
poliovirus [@PMID-17223130] which do contain nucleotides to map residues that
are expected to interact with nucleotides. The mapping is possible since despite
remote homology, the polymerase architecture is largely conserved and thus
allows for structural superpositions.

![Differences to SARS-CoV (red) vs predicted ATP binding site (blue) on RNA polymerase \label{RdRp_combined}](./RdRp-combined.png)

## Annotations extracted from literature

We generated annotations based on the most recent literature on the SARS-CoV-2 spike protein (P0DTC2) receptor binding domain (RBD). This included putative adaptive mutations from SARS-CoV to SARS-CoV-2 and highly conserved contact residues in the RBD [@Wan2020; @Tai2020]. Residues related to ACE2 binding [@Shang2020] and CR3022 antibody interactions were also highlighted [@Yuan2020]. Mapping annotations extracted from recent literature onto the same protein structure can allow for comparison and synthesis of information across research groups and the formation of new hypotheses. These annotations can also be used to add detail to those extracted from variations within sequence data.

# Conclusion

TODO:
- Some conclusion summarizing the results

# Future work

TODO:
- Next steps?

# GitHub repository for data and scripts

* All annotations and scripts are stored here: https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures

# Acknowledgements
This work was done within the [Covid-19 Virtual BioHackathon 2020](https://github.com/virtual-biohackathons/covid-19-bh20).
We thank the organizers for the opportunity and the support.

# References

