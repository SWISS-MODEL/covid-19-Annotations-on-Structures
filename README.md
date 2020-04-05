# Mapping sequence data onto structures
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-5-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

This repository collects contributions related to the ["Annotations on Structures" topic](https://github.com/virtual-biohackathons/covid-19-bh20/wiki/Annotations-on-Structures) in the [COVID-19 Biohackathon April 5-11 2020](https://github.com/virtual-biohackathons/covid-19-bh20).

The context is [SWISS-MODEL's](https://swissmodel.expasy.org) involvement in an [EU project to combat COVID-19](https://www.sib.swiss/about-sib/news/10659). To accelerate our plan to map relevant annotations onto those structures, we collect tools/platforms which can automatically generate such annotations based on the latest data.

We mainly hope to receive two types of contributions:
1. Find/generate relevant sequence data (see [issues list](https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures/issues) for inspirational ideas) to be displayed on structures (see [section on SWISS-MODEL's annotation system](#swiss-model-annotation-system)). This should be scripted to enable automated fetching of the latest data.
2. Write reusable scripts to map the sequence data onto the frame of reference of proteins (this might need translation from position on genome data to position on proteins of SARS-CoV-2 as listed [here](https://swissmodel.expasy.org/repository/species/2697049)). These scripts are expected to be useful for the scripts in point 1.

Additional topics of interest:
- For visualization experts: alternative ways to visualize the protein structures.
- For RDF/JSON-LD experts: define an RDF ontology and map our json-data ([example](https://swissmodel.expasy.org/repository/uniprot/P59594.json)) to RDF to be used in other knowledge graph efforts. Some efforts exist from [PDBj](https://pdbj.org/help/rdf) to map structures to RDF but they focus on experimental meta data while we consider structural coverage of the proteins more relevant. Probably [SIFTS mappings](https://pdbj.org/news/20160629) are the better starting point here. With a minimal "@context" section referring to UniProt we might also be able to turn our existing json to valid json-ld.
- For protein modelling experts: custom modeling of proteins of interest (e.g. using careful expert-curated target-template alignments or combination of templates)

## Preferred technologies

- Programming languages used within SWISS-MODEL: Python (3.6), C++
- Dealing with protein structure and sequence data: [OpenStructure](https://openstructure.org/)

## Guidelines for contributions

Follow the biohackathon's [code of conduct](https://github.com/virtual-biohackathons/covid-19-bh20/blob/master/CODE_OF_CONDUCT.md) and this project's [contributions guidelines](CONTRIBUTING.md).

## SWISS-MODEL annotation system

**NOTE: this is work-in-progress and subject to change.**

The beta-server of SWISS-MODEL is used to allow users to upload annotations: https://beta.swissmodel.expasy.org/repository/covid_annotation_upload

Both the user annotations and the display of the viral polyprotein ([R1AB_SARS2](https://beta.swissmodel.expasy.org/repository/uniprot/P0DTD1)) are still work-in-progress and may have bugs. If you find problems with those prototype SWISS-MODEL features, please add issues to this github project and we will try to address them as soon as possible.

The annotation format is a plain-text format:
- One line per annotation
- Each annotation will consist of 5 or 6 comma- or tab-separated values:
  1. ID (UniProtKB AC or MD5 checksum of the sequence)
  2. Start position (1-based)
  3. End position
  4. Color value
  5. Reference (optional)
  6. Annotation comment
- Example:
  ```
  P0DTD1	3400	3450	#FF00FF	https://swissmodel.expasy.org/repository/	My Awesome Annotation
  P0DTC2	230	330	#FFA500	A text reference	One more!
  ```
- UniProtKB ACs with links can be found in [UniProtKB](https://covid-19.uniprot.org/)
  - Our [SARS-CoV-2 page](https://swissmodel.expasy.org/repository/species/2697049) shows mapping to mature proteins and the correspondance to RefSeq and GenBank.
  - For cleaved proteins, use the parent protein. For instance an annotation on nsp3 (Non-structural protein 3) must be reported on P0DTD1 (the "parent" protein) with an offset of 818 (as nsp3 start on position 819 of P0DTD1).

Also we are actively working on extending the structural coverage of the SARS-CoV-2 proteome by using protein predictions from colleagues participating in CASP.

## Context

Protein structure predictions of SARS-CoV-2 have already proven useful to several research projects. To list a few examples which used our models:
- [A potential role for integrins in host cell entry by SARS-CoV-2, Antiviral Research](https://doi.org/10.1016/j.antiviral.2020.104759)
- [Targeting Novel Coronavirus 2019: A Systematic Drug Repurposing Approach to Identify Promising Inhibitors Against 3C-like Proteinase and 2'-O-Ribose Methyltransferase](https://dx.doi.org/10.26434/chemrxiv.11888730.v1)
- [Genomic characterisation and epidemiology of 2019 novel coronavirus: implications for virus origins and receptor binding, The Lancet](https://dx.doi.org/10.1016/S0140-6736(20)30251-8)
- [Insilico Medicine publishes molecular structures for the key protein target of 2019-nCoV](https://insilico.com/ncov-sprint)
- [Targeting 2019-nCoV: GHDDI Info Sharing Portal](https://ghddi-ailab.github.io/Targeting2019-nCoV/)

## Contributors ✨

Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/gtauriello"><img src="https://avatars3.githubusercontent.com/u/25968022?v=4" width="100px;" alt=""/><br /><sub><b>Gerardo Tauriello</b></sub></a><br /><a href="#projectManagement-gtauriello" title="Project Management">📆</a></td>
    <td align="center"><a href="https://github.com/xrobin"><img src="https://avatars2.githubusercontent.com/u/1047170?v=4" width="100px;" alt=""/><br /><sub><b>Xavier Robin</b></sub></a><br /><a href="#tool-xrobin" title="Tools">🔧</a> <a href="https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures/commits?author=xrobin" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/bienchen"><img src="https://avatars0.githubusercontent.com/u/69343?v=4" width="100px;" alt=""/><br /><sub><b>bienchen</b></sub></a><br /><a href="#tool-bienchen" title="Tools">🔧</a></td>
    <td align="center"><a href="https://github.com/awaterho"><img src="https://avatars2.githubusercontent.com/u/40768716?v=4" width="100px;" alt=""/><br /><sub><b>Andrew W</b></sub></a><br /><a href="#tool-awaterho" title="Tools">🔧</a> <a href="#design-awaterho" title="Design">🎨</a></td>
    <td align="center"><a href="https://github.com/schdaude"><img src="https://avatars3.githubusercontent.com/u/4851123?v=4" width="100px;" alt=""/><br /><sub><b>schdaude</b></sub></a><br /><a href="#tool-schdaude" title="Tools">🔧</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
