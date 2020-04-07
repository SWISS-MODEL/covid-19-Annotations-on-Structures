# Mapping sequence data onto structures
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-12-orange.svg?style=flat-square)](#contributors-)
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
- Dealing with protein structure and sequence data: [OpenStructure](https://openstructure.org/) (example in [wiki here](https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures/wiki/Annotation-example-with-the-OpenStructure-Computational-Structural-Biology-Framework))

## Guidelines for contributions

Follow the biohackathon's [code of conduct](https://github.com/virtual-biohackathons/covid-19-bh20/blob/master/CODE_OF_CONDUCT.md) and this project's [contributions guidelines](CONTRIBUTING.md).

## SWISS-MODEL annotation system

**NOTE: this is work-in-progress and subject to change.**

The beta-server of SWISS-MODEL is used to allow users to upload annotations: https://beta.swissmodel.expasy.org/repository/covid_annotation_upload (a list of projects for registered users can be found [here](https://beta.swissmodel.expasy.org/repository/covid_annotation_projects)).

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
- The Annotation class available in utils facilitates creation of new annotations:
  ```python
  from utils.sm_annotations import Annotation

  # generate example annotations
  annotation = Annotation()

  # Annotation of residue range with color red provided as RGB
  annotation.add("P0DTD1", (10, 20), (1.0, 0.0, 0.0), "red anno")

  # Again, annotating a range but this time we're adding a reference
  # and provide the color blue as hex
  annotation.add("P0DTD1", (21, 30), "#0000FF", "blue anno", 
                 reference = "https://swissmodel.expasy.org/")

  # Outputs plain text which is accepted on the covid annotation upload 
  print(annotation)

  # Or directly do a post request (defaults to SWISS-MODEL beta)
  print("visit", annotation.post(), "to see awesome things")
  ```
  The last line directly creates a new annotation project and prints its url. 
  An example can be viewed [here](https://beta.swissmodel.expasy.org/repository/covid_annotation_project/dUHPPN)

- UniProtKB ACs with links can be found in [UniProtKB](https://covid-19.uniprot.org/)
  - Our [SARS-CoV-2 page](https://swissmodel.expasy.org/repository/species/2697049) shows mapping to mature proteins and the correspondence to RefSeq and GenBank.
  - We also have a [list of all SARS-CoV-2 proteins](https://beta.swissmodel.expasy.org/repository/species/2697049/list) that shows an overview of the ACs and their structural coverage.
  - For cleaved proteins, use the parent protein. For instance an annotation on nsp3 (Non-structural protein 3) must be reported on P0DTD1 (the "parent" protein) with an offset of 818 (as nsp3 start on position 819 of P0DTD1).
  - ViralZone has a well described overview of the proteome [here](https://viralzone.expasy.org/8996).
  - We propose to ignore the shorter polyprotein (P0DTC1, R1A_SARS2) as it's cleaved into the same mature proteins as the longer one (P0DTD1, R1AB_SARS2) with the exception of a very short peptide (Non-structural protein 11 (nsp11), [YP_009725312.1](https://www.ncbi.nlm.nih.gov/protein/YP_009725312.1)).
  - Two proteins of unknown function ([P0DTD2](https://covid-19.uniprot.org/uniprotkb/P0DTD2) and [P0DTD3](https://covid-19.uniprot.org/uniprotkb/P0DTD3)) are missing from our [SARS-CoV-2 page](https://swissmodel.expasy.org/repository/species/2697049) but can safely be used to map annotations and we will provide structures if possible.
  - Additionally to the SARS-CoV-2 proteins, it also makes sense to map annotations for [Q9BYF1](https://covid-19.uniprot.org/uniprotkb/Q9BYF1) (ACE2_HUMAN). So far this is the only virus-host-interaction for which we have structural information. More interactions have been proposed (e.g. [here](https://viralzone.expasy.org/9077)) but we don't have structures for them (yet).

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
    <td align="center"><a href="https://github.com/schdaude"><img src="https://avatars3.githubusercontent.com/u/4851123?v=4" width="100px;" alt=""/><br /><sub><b>schdaude</b></sub></a><br /><a href="#tool-schdaude" title="Tools">🔧</a> <a href="https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures/commits?author=schdaude" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/BarbaraTerlouw"><img src="https://avatars0.githubusercontent.com/u/47810869?v=4" width="100px;" alt=""/><br /><sub><b>BarbaraTerlouw</b></sub></a><br /><a href="#ideas-BarbaraTerlouw" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/vprobon"><img src="https://avatars1.githubusercontent.com/u/49338525?v=4" width="100px;" alt=""/><br /><sub><b>Vasilis J Promponas</b></sub></a><br /><a href="#ideas-vprobon" title="Ideas, Planning, & Feedback">🤔</a></td>
  </tr>
  <tr>
    <td align="center"><a href="http://biohackathons.github.io"><img src="https://avatars0.githubusercontent.com/u/5738421?v=4" width="100px;" alt=""/><br /><sub><b>Ben Busby</b></sub></a><br /><a href="#ideas-DCGenomics" title="Ideas, Planning, & Feedback">🤔</a> <a href="#content-DCGenomics" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/lnblum"><img src="https://avatars2.githubusercontent.com/u/51452159?v=4" width="100px;" alt=""/><br /><sub><b>Laura Blum</b></sub></a><br /><a href="#content-lnblum" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/tomasMasson"><img src="https://avatars0.githubusercontent.com/u/59352285?v=4" width="100px;" alt=""/><br /><sub><b>tomasMasson</b></sub></a><br /><a href="#content-tomasMasson" title="Content">🖋</a> <a href="https://github.com/SWISS-MODEL/covid-19-Annotations-on-Structures/commits?author=tomasMasson" title="Code">💻</a></td>
    <td align="center"><a href="http://www.linkedin.com/in/didier-barradas-bautista"><img src="https://avatars3.githubusercontent.com/u/17081199?v=4" width="100px;" alt=""/><br /><sub><b>Didier Barradas Bautista</b></sub></a><br /><a href="#content-D-Barradas" title="Content">🖋</a></td>
    <td align="center"><a href="https://github.com/bmeldal"><img src="https://avatars2.githubusercontent.com/u/10517124?v=4" width="100px;" alt=""/><br /><sub><b>Birgit Meldal</b></sub></a><br /><a href="#ideas-bmeldal" title="Ideas, Planning, & Feedback">🤔</a> <a href="#content-bmeldal" title="Content">🖋</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
