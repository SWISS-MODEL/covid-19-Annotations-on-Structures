import csv
from argparse import ArgumentParser
from utils.uniprot import seq_from_ac
from utils.sm_annotations import Annotation


def _parse_args():
    parser = ArgumentParser(
        description="Loads csv generated by get_nextstrain_data.py and "
        "collects annotations that can be displayed in the SWISS-MODEL "
        "annotation system."
    )
    parser.add_argument(
        "--csv",
        required=True,
        type=str,
        dest="csv",
        help="CSV file generated by get_nextstrain_data.py",
    )
    parser.add_argument(
        "--post",
        dest="post",
        help="If set, the annotations "
        "are directly uploaded to the SWISS-MODEL annotation system and an "
        "accessible URL is printed to stdout. By default, the annotations are "
        "printed to stdout to be copy pasted in the annotation upload form in "
        "SWISS-MODEL",
        action="store_true",
    )
    return parser.parse_args()


def _parse_mutations(csv_file):
    mutations = list()
    with open(csv_file, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            mutations.append((row["Protein"], row["Mutation"]))
    return mutations


def _annotate(mutations):
    """
    The nextstrain data is relative to a phylogenetic tree that they build.
    So a mutation at location x can either be a mutation from reference 
    to something else, a mutation from something else back to the reference
    or even something that has nothing to do with the reference at all.
    So the variations we're reporting are all observed amino acid one letter 
    codes (olc) at a certain location => the unique set of olcs containing the 
    the reference and all corresponding olcs reported in the nextstrain 
    mutations.
    """

    # lets fetch all reference sequences
    acs = [m[0] for m in mutations]
    sequences = dict()
    for ac in set(acs):
        sequences[ac] = seq_from_ac(ac)

    # key: (uniprot_ac, rnum), value: [olc1, olc2, ..., olcx]
    variations = dict()
    for m in mutations:
        uniprot_ac = m[0]
        rnum = int(m[1][1 : len(m[1]) - 1])
        key = (uniprot_ac, rnum)
        if key not in variations:
            variations[key] = set()
        variations[key].add(sequences[uniprot_ac][rnum - 1])
        variations[key].add(m[1][0])  # first letter -> mutation from
        variations[key].add(m[1][-1])  # last letter -> mutation to

    # Let's not directly add it to the annotation object but rather do an
    # intermediate list so we can sort it by uniprot accession code / residue
    # number
    annotation_data = list()
    for var_key in variations.keys():
        annotation_data.append(
            [var_key[0], var_key[1], "r", ",".join(variations[var_key])]
        )
    annotation = Annotation()
    for d in sorted(annotation_data):
        annotation.add(d[0], d[1], d[2], d[3])

    return annotation


def main():
    args = _parse_args()
    mutations = _parse_mutations(args.csv)
    annotation = _annotate(mutations)
    if args.post:
        print(annotation.post())
    else:
        print(annotation)


if __name__ == "__main__":
    main()
