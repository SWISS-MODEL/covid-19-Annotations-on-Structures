# Mapping sequence data onto structures

The [SWISS-MODEL](https://swissmodel.expasy.org) team is currently involved in an [EU project to combat COVID-19](https://www.sib.swiss/about-sib/news/10659).

There we are providing protein structures as starting point for further analysis (see [here](https://swissmodel.expasy.org/repository/species/2697049)). Given the current outbreak, we would now like to accelerate our plan to map relevant annotations onto those structures. Hence, we are very much interested in tools/platforms which can automatically generate such annotations based on the latest data.

Within this hackathon, the goals are to:
1. Find/generate relevant sequence data (e.g. variations from different strains) to be mapped onto structures (this should be automatized to always fetch the latest data)
2. Write scripts to map the sequence data onto the frame of reference of proteins (this might need translation from position on genome data to position on proteins of SARS-CoV-2 as listed [here](https://swissmodel.expasy.org/repository/species/2697049))

We expect that point 1 can have many parallel developments and that point 2 can be reused for different sequence data sources.

We will provide an interface to display annotations in a similar fashion as done in our SWISS-MODEL repository (as in [this example](https://swissmodel.expasy.org/repository/uniprot/B8XC04)). That being said, if anyone wants to work on alternative ways to visualize the protein structures, we are happy to provide the structures. Also we will be extending the structural coverage of the SARS-CoV-2 proteome by using protein predictions from colleagues participating in CASP.

## Preferred technologies

- Programming languages used within SWISS-MODEL: Python (3.6), C++
- Dealing with protein structure and sequence data: [OpenStructure](https://openstructure.org/)

## Guidelines for contributions

Follow the biohackathon's [code of conduct](https://github.com/virtual-biohackathons/covid-19-bh20/blob/master/CODE_OF_CONDUCT.md) and this project's [contributions guidelines](CONTRIBUTING.md).

## SWISS-MODEL annotation system

**NOTE: this is work-in-progress and subject to change.**

The beta-server of SWISS-MODEL will be used to allow user-annotations to be uploaded (more details to follow once this is enabled).

The annotation format is a plain-text format:
- One line per annotation
- Each annotation will consist of 5 or 6 space-, comma- or tab-separated values:
  1. ID (UniProtKB AC or MD5 checksum of the sequence)
  2. Start position (1-based)
  3. End position
  4. Color value
  5. Reference (optional)
  6. Annotation comment
- Example:
  ```
  P0DTD1	3400	3450	#FF00FF	https://swissmodel.expasy.org/repository/	My Awesome Annotation
  P0DTC2	230	330	#FFA500	One more!
  ```
- UniProtKB ACs with links can be found in our [SARS-CoV-2 page](https://swissmodel.expasy.org/repository/species/2697049)

## Context

Protein structure predictions of SARS-CoV-2 have already proven useful to several research projects. To list a few examples which used our models:
- [A potential role for integrins in host cell entry by SARS-CoV-2, Antiviral Research](https://doi.org/10.1016/j.antiviral.2020.104759)
- [Targeting Novel Coronavirus 2019: A Systematic Drug Repurposing Approach to Identify Promising Inhibitors Against 3C-like Proteinase and 2'-O-Ribose Methyltransferase](https://dx.doi.org/10.26434/chemrxiv.11888730.v1)
- [Genomic characterisation and epidemiology of 2019 novel coronavirus: implications for virus origins and receptor binding, The Lancet](https://dx.doi.org/10.1016/S0140-6736(20)30251-8)
- [Insilico Medicine publishes molecular structures for the key protein target of 2019-nCoV](https://insilico.com/ncov-sprint)
- [Targeting 2019-nCoV: GHDDI Info Sharing Portal](https://ghddi-ailab.github.io/Targeting2019-nCoV/)
