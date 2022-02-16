# Fetch variation data from Nextstrain API

First draft script to retrieve and parse the mutations reported in raw data from [Nextstrain](http://data.nextstrain.org/ncov.json). 

I include all the variation data (amino acid and nucleotide) and all the metadata (we could decide later which fields are necessary).

Annotations describing all distinct amino acid one letter codes or all mutation 
events occuring at a certain location can be fetched from the generated 
nextstrain_data.csv and prepared for the SWISS-MODEL annotation system by 
the annotate_variation.py script.

To run it, you require the the root directory of the github 
repository in your PYTHONPATH. On Linux/Mac execute the following:

```
export PYTHONPATH=<path_to_covid_repo>:$PYTHONPATH
```

Usage and available options of the script can be displayed by executing

```
python annotate_variation.py --help
```

- Example annotation of all distinct amino acids is available [here](https://beta.swissmodel.expasy.org/repository/covid_annotation_project/LPRPYc)
- Example annotation of all observed mutation events is available [here](https://beta.swissmodel.expasy.org/repository/covid_annotation_project/vtvJJN)

TODO:

- Give format to the output (done).
- Select the fields of interest.
- Replace protein names with UniProtKB AC (done).

