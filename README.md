# Mapping sequence data onto structures
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-2-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->

This repository collects contributions related to the ["Annotations on Structures" topic](https://github.com/virtual-biohackathons/covid-19-bh20/wiki/Annotations-on-Structures) in the [COVID-19 Biohackathon April 5-11 2020](https://github.com/virtual-biohackathons/covid-19-bh20).

The context is [SWISS-MODEL's](https://swissmodel.expasy.org) involvement in an [EU project to combat COVID-19](https://www.sib.swiss/about-sib/news/10659). To accelerate our plan to map relevant annotations onto those structures, we collect tools/platforms which can automatically generate such annotations based on the latest data.

We mainly hope to receive two types of contributions:
1. Find/generate relevant sequence data (see [ideas below](#ideas-for-annotations)) to be displayed on structures (see [section on SWISS-MODEL's annotation system](#swiss-model-annotation-system)). This should be scripted to enable automated fetching of the latest data.
2. Write reusable scripts to map the sequence data onto the frame of reference of proteins (this might need translation from position on genome data to position on proteins of SARS-CoV-2 as listed [here](https://swissmodel.expasy.org/repository/species/2697049)). These scripts are expected to be useful for the scripts in point 1.

## Ideas for annotations

For inspiration, here are some annotation ideas with the expected work to be done:
- Include variations from processed data in [nextstrain](https://nextstrain.org/ncov). Requires:
  - parse [json data](https://data.nextstrain.org/ncov.json) following their [dev docs](https://github.com/nextstrain/ncov/blob/master/DEV_DOCS.md)
  - map variations onto UniProtKB ACs used in [SWISS-MODEL](https://swissmodel.expasy.org/repository/species/2697049)
  - define colors and annotation texts variations
  - test using [SWISS-MODEL's annotation system](#swiss-model-annotation-system)
  - properly acknowledge source of data (see also "Data" section in nextstrain's [README](https://github.com/nextstrain/ncov/blob/master/README.md))
  - followups: add possibility to filter results (e.g. only from country X or certain confidence), process into entropies, ...
- Process glycosylation sites and map expected changes of accessibility onto the structures (see [here](https://twitter.com/Olivercgrant/status/1243576788514725888), [here](https://twitter.com/rommieamaro/status/1241810976866840577?s=11), [here](https://twitter.com/ElisaTelisa/status/1244174688437374978) and [here](https://www.biorxiv.org/content/10.1101/2020.03.28.013276v1.full.pdf) for work on this).
- Go beyond the UniProt-annotations which we already display. E.g. for PL-PRO of SARS-CoV we found annotations in the literature (see Fig. 1 of [Báez-Santos et al 2015](https://doi.org/10.1016/j.antiviral.2014.12.015)) which go way beyond what we find in the [UniProt-annotations](https://covid-19.uniprot.org/).

This list is by no means exhaustive and we can imagine that a lot more useful annotations can be found by parsing the ever-increasing literature on SARS-CoV-2.

## Preferred technologies

- Programming languages used within SWISS-MODEL: Python (3.6), C++
- Dealing with protein structure and sequence data: [OpenStructure](https://openstructure.org/)

## Guidelines for contributions

Follow the biohackathon's [code of conduct](https://github.com/virtual-biohackathons/covid-19-bh20/blob/master/CODE_OF_CONDUCT.md) and this project's [contributions guidelines](CONTRIBUTING.md).

## SWISS-MODEL annotation system

**NOTE: this is work-in-progress and subject to change.**

The beta-server of SWISS-MODEL will be used to allow user-annotations to be uploaded (more details to follow once this is enabled). The annotations will be displayed in a similar fashion as done in our SWISS-MODEL repository (as in [this example](https://swissmodel.expasy.org/repository/uniprot/B8XC04)).

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

If anyone wants to work on alternative ways to visualize the protein structures, we are happy to provide the structures. Also we are actively working on extending the structural coverage of the SARS-CoV-2 proteome by using protein predictions from colleagues participating in CASP.

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
    <td align="center"><a href="https://github.com/xrobin"><img src="https://avatars2.githubusercontent.com/u/1047170?v=4" width="100px;" alt=""/><br /><sub><b>Xavier Robin</b></sub></a><br /><a href="#tool-xrobin" title="Tools">🔧</a></td>
  </tr>
</table>

<!-- markdownlint-enable -->
<!-- prettier-ignore-end -->
<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!