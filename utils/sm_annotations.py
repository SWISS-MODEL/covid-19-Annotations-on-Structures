import matplotlib as mpl
import re
import requests
from . import uniprot


class Annotation:
    """
    Helper class to programmatically define annotations and format according
    to SWISS-MODEL requirements

    usage example:

    from sm_annotations import Annotation

    # generate example annotations
    annotation = Annotation()

    # Annotation of residue range with color red provided as RGB
    annotation.add("P0DTD1", (10, 20), (1.0, 0.0, 0.0), "red anno")

    # Again, annotating a range but this time we're adding a reference
    # and provide the color blue as hex
    annotation.add("P0DTD1", (21, 30), "#0000FF", "blue anno", 
                   reference = "https://swissmodel.expasy.org/")

    # Single residue annotation in green
    annotation.add("P0DTD1", 35, "#00FF00", "green anno")

    # print string which is accepted on 
    # https://swissmodel.expasy.org/repository/user_annotation_upload
    print(annotation)

    # or directly do a post request (defaults to SWISS-MODEL)
    print("Visit the following url to see awesome things:")
    print(annotation.post(title="awesome things"))
    """

    def __init__(self):
        self.uniprot_acs = list()
        self.rnum_ranges = list()
        self.colors = list()
        self.annotations = list()
        self.references = list()

    def add(self, uniprot_ac, rnum, color, annotation, reference=None):
        """
        Check for valid data and add new annotation

        :param uniprot_ac:   Valid UniprotAC as string to which the annotation 
                             referes to
        :param rnum:         Specify location of annotation. 
                             Can either be an integer describing a single 
                             residue annotation or a tuple/list with two 
                             elements specifying a full range. Numbers refer
                             to the specified Uniprot reference sequence with
                             numbering starting at 1.
        :param color:        Color of annotation. Can be any "color-like" 
                             expression accepted by matplotlib
                             (https://matplotlib.org/3.2.1/api/colors_api.html)
        :param annotation:   The actual annotation as string
        :param reference:    Optional reference as string, e.g. URL or free text. 
        """

        # check input
        if not uniprot.valid_uniprot_ac_pattern(uniprot_ac):
            raise ValueError("uniprot_ac is invalid")

        if isinstance(rnum, int):
            if rnum < 1:
                raise ValueError("Expect rnum >= 1")
        elif isinstance(rnum, tuple) or isinstance(rnum, list):
            if len(rnum) != 2:
                raise ValueError(
                    "Expect rnum to be integer or tuple/list with two elements"
                )
            if rnum[0] < 1 or rnum[1] < 1:
                raise ValueError("Expect rnum >= 1 for all rnum")

        if not mpl.colors.is_color_like(color):
            raise ValueError(
                "Only accept color formats specified in https://matplotlib.org/3.2.1/api/colors_api.html"
            )

        if not isinstance(annotation, str):
            raise ValueError("Expect annotation to be str")

        if reference:
            if not isinstance(reference, str):
                raise ValueError("Expect reference to be None or str")

        # all input is valid, add annotation
        self.uniprot_acs.append(uniprot_ac)
        if isinstance(rnum, int):
            self.rnum_ranges.append((rnum, rnum))
        else:
            self.rnum_ranges.append(rnum)
        self.colors.append(color)
        self.annotations.append(annotation)
        self.references.append(reference)

    def post(
        self,
        url="https://swissmodel.expasy.org/repository/user_annotation_upload",
        title=None,
        email=None,
    ):
        """
        Performs post request and returns the url at which the annotations can
        be accessed. Usage example:

        print("visit", annotation.post(), "to see awesome things")

        :param url:     URL of annotation upload form, defaults to SWISS-MODEL.
        :param title:   Filled in project title field if given
        :param email:   Filled in email field if given
        """
        data = {"annotation_data": self.__str__()}
        if title:
            data["title"] = title
        if email:
            data["email"] = email
        res = requests.post(url, data=data)
        # Ensure we didn't get an error
        res.raise_for_status()
        # Ensure we were redirected
        if res.url == url:
            raise ValueError("The submission of annotations failed")
        return res.url

    def __str__(self):
        """
        Processes annotation and return string which is accepted on
        https://swissmodel.expasy.org/repository/user_annotation_upload
        """
        n = len(self.uniprot_acs)
        formatted_annos = [self._format_anno(idx) for idx in range(n)]
        return "\n".join(formatted_annos)

    def __len__(self):
        return len(self.annotations)

    def _format_anno(self, idx):
        if idx < 0 or idx >= len(self.uniprot_acs):
            raise ValueError("Invalid idx in line formatting")
        data = list()
        data.append(self.uniprot_acs[idx])
        data.append(str(self.rnum_ranges[idx][0]))
        data.append(str(self.rnum_ranges[idx][1]))
        data.append(mpl.colors.to_hex(self.colors[idx]))
        if self.references[idx]:
            data.append(self.references[idx])
        data.append(self.annotations[idx])
        return "\t".join(data)
