from utils.sm_annotations import Annotation
from utils import uniprot


def _align(s_one, s_two):
    # for a hacky first implementation we rely on the OpenStructure package
    # would be cool to avoid that dependency for this simple example annotation
    from ost import seq

    ost_s_one = seq.CreateSequence("one", s_one)
    ost_s_two = seq.CreateSequence("two", s_two)
    aln = seq.alg.GlobalAlign(ost_s_one, ost_s_two, seq.alg.BLOSUM62)[0]
    return (str(aln.GetSequence(0)), str(aln.GetSequence(1)))


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
    aligned_s_target, aligned_s_reference = _align(s_target, s_reference)

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
