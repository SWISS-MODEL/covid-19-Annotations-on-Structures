from utils.sm_annotations import Annotation
from utils import uniprot
from needleman_wunsch import Align


def get_annotation(uniprot_ac_target, uniprot_ac_reference, nonconserved_color="r"):
    """
    Annotate mutations of a target sequence relative to a reference 
    sequence. The alignment of the sequences is performed using a pairwise
    Needleman-Wunsch algorithm.
    
    :param uniprot_ac_target: Sequence for which you want the annotations
    :param uniprot_ac_reference: Reference sequence from which changes are
                                 annotated
    :param nonconserved_color: Color assigned to annotations of non-conserved
                               amino acids
    """
    s_target = uniprot.seq_from_ac(uniprot_ac_target)
    s_reference = uniprot.seq_from_ac(uniprot_ac_reference)
    aligned_s_target, aligned_s_reference = Align(s_target, s_reference)

    # this is a sanity check if the align function did not alter the underlying
    # sequenes
    assert s_target == aligned_s_target.replace("-", "")
    assert s_reference == aligned_s_reference.replace("-", "")
    assert len(aligned_s_target) == len(aligned_s_reference)

    target_rnum = 0
    annotation = Annotation()
    for i in range(len(aligned_s_target)):
        if aligned_s_target[i] != "-":
            target_rnum += 1
        if aligned_s_target[i] != "-" and aligned_s_reference[i] != "-":
            if aligned_s_target[i] != aligned_s_reference[i]:
                annotation.add(
                    uniprot_ac_target,
                    target_rnum,
                    nonconserved_color,
                    "%s->%s" % (aligned_s_reference[i], aligned_s_target[i]),
                )

    return annotation
