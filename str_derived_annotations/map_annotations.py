from utils.sm_annotations import Annotation
import matplotlib as mpl
from matplotlib import cm
import pandas as pnd
import numpy as np

def numbers_to_colors(numbers, cmap="jet", log=False):
    if log:
        numbers = np.log1p(numbers)
    norm = mpl.colors.Normalize(vmin=np.min(numbers), vmax=np.max(numbers))
    colormap = cm.get_cmap(cmap)
    return [colormap(norm(n))[:3] for n in numbers]

def example_annotations_ensemble():
    main_protease_id = "P0DTD1"
    reference_pdb_id = "6y2g"
    offset = 3264
    df = pnd.read_csv("annotations/6y2g_AAA_reference.txt", sep="\t")
   
    
    mapper = {}
    for c in ["RMSDs per residue", "PCA fluctuations"]:
        mapper[c] = numbers_to_colors(df[c].values)
    for c in mapper:
        annotation = Annotation()
        for index, row in df.iterrows():
            annotation.add(main_protease_id, row["Residue Index"] + offset, tuple(mapper[c][index]), c, reference=reference_pdb_id)
        print(annotation.post(title=f"{c}_{reference_pdb_id}"))
        
        
def example_annotations_anm():
    main_protease_id = "P0DTD1"
    reference_pdb_id = "6y2g"
    offset = 3264
    df = pnd.read_csv("annotations/6y2g_AAA.txt", sep="\t")
    mapper = {}
    for c in df.columns:
        if c == "Index":
            continue
        mapper[c] = numbers_to_colors(df[c].values)
    for c in mapper:
        annotation = Annotation()
        for index, row in df.iterrows():
            annotation.add(main_protease_id, int(row["Index"]) + offset, tuple(mapper[c][index]), c, reference=reference_pdb_id)
        print(annotation.post(title=f"{c}_{reference_pdb_id}"))

if __name__ == "__main__":
    example_annotations_anm()