import numpy as np

"""
Code that I had lying around... I hope it works as desired. I'm fully aware 
that it cannot keep up with more efficient Needleman-Wunsch implementations.
However, we have dependency free code! yeah!
"""

class SubstitutionMatrix:
    def __init__(self, data_file="blosum62.txt"):

        data = open(data_file, "r").readlines()

        if len(data) != 21:
            err = "Expect exaclty 21 lines in matrix file!"
            err += " One line defining the amino acid one letter codes and "
            err += "20 lines defining the substitution matrix"
            raise RuntimeError(err)

        self.one_letter_codes = data[0].strip()

        if len(self.one_letter_codes) != 20:
            err = "First line must contain exactly 20 letters representing "
            err += "the amino acids!"
            raise RuntimeError(err)

        self.matrix = np.zeros((20, 20))

        for i in range(20):
            data_line = data[1 + i].split()
            if len(data_line) != 20:
                raise RuntimeError("Each data line must contain exactly 20 elements!")
            for j in range(20):
                self.matrix[(i, j)] = int(data_line[j])

    def GetScore(self, aa_one, aa_two):
        """
    Substitution score given two one letter codes

    :param aa_one:      one letter code of first amino acid
    :param aa_two:      one letter code of second amino acid
    """

        idx_one = self.one_letter_codes.find(aa_one)
        idx_two = self.one_letter_codes.find(aa_two)

        if idx_one == -1:
            raise RuntimeError("aa_one must be one of: " + self.one_letter_codes)

        if idx_two == -1:
            raise RuntimeError("aa_two must be one of: " + self.one_letter_codes)

        return self.matrix[(idx_one, idx_two)]


def Align(s1, s2, gap_penalty=-8, subst_matrix=None):
    """
    Aligns two raw strings using a Needleman-Wunsch algorithm and returns a
    tuple containing the aligned input. '-' represent gaps.

    :param s1:       String representing the first sequence
    :param s2:       String representing the second sequence
    :param gap_penalty: Penalty value for opening/extending a gap
    :param subst_matrix: SubstitutionMatrix object for scoring, 
                         defaults to BLOSUM62 parametrization.
    """

    if subst_matrix is None:
        subst_matrix = SubstitutionMatrix()

    # setup scoring and backtracking matrix
    n_rows = len(s1) + 1  # s1 goes from top to bottom
    n_cols = len(s2) + 1  # s2 goes from left to right
    scoring_matrix = np.zeros((n_rows, n_cols))

    # backtrack encoding:
    # 1 => aligned
    # 2 => deletion in s1 (insertion in s2 respectively)
    # 3 => deletion in s2 (insertion in s1 respectively)
    backtrack_matrix = np.zeros((n_rows, n_cols))

    # the first row and the first column can already be prefilled
    for i in range(1, n_rows):
        scoring_matrix[(i, 0)] = i * gap_penalty
        backtrack_matrix[(i, 0)] = 3
    for i in range(1, n_cols):
        scoring_matrix[(0, i)] = i * gap_penalty
        backtrack_matrix[(0, i)] = 2

    # fill scoring and backtracking matrices
    #
    # Every position in the scoring matrix represents the best local solution given the
    # steps I can take
    # - Align two residues => The value of the scoring matrix to the upper left plus the
    #                         according score from the scoring matrix
    # - Apply deletion in s1 => The value from the scoring matrix to the left plus
    #                           the penalty for a gap
    # - Apply deletion in s2 => The value from the scoring matrix above plus
    #                           the penalty for a gap
    #
    # The backtracking matrix stores the path I took, e.g. 1 for the first option
    for r_idx in range(1, n_rows):
        for c_idx in range(1, n_cols):
            aligned_score = scoring_matrix[
                (r_idx - 1, c_idx - 1)
            ] + subst_matrix.GetScore(s1[r_idx - 1], s2[c_idx - 1])
            s1_deletion_score = scoring_matrix[(r_idx, c_idx - 1)] + gap_penalty
            s2_deletion_score = scoring_matrix[(r_idx - 1, c_idx)] + gap_penalty
            scoring_matrix[(r_idx, c_idx)] = max(
                aligned_score, s1_deletion_score, s2_deletion_score
            )

            if aligned_score > s1_deletion_score and aligned_score > s2_deletion_score:
                backtrack_matrix[(r_idx, c_idx)] = 1
            elif s1_deletion_score > s2_deletion_score:
                backtrack_matrix[(r_idx, c_idx)] = 2
            else:
                backtrack_matrix[(r_idx, c_idx)] = 3

    # perform backtracking to get final alignment
    # In principle we start at the bottom right of the backtracking matrix and
    # work our way through the matrix until we hit the upper left
    r_idx = n_rows - 1
    c_idx = n_cols - 1
    path = list()

    # if we hit zero, we're at the upper left corner and can stop
    while backtrack_matrix[(r_idx, c_idx)] != 0:
        path.append(backtrack_matrix[(r_idx, c_idx)])
        if backtrack_matrix[(r_idx, c_idx)] == 1:
            r_idx -= 1
            c_idx -= 1
        elif backtrack_matrix[(r_idx, c_idx)] == 2:
            c_idx -= 1
        else:
            r_idx -= 1

    # backtracking comes from the back, so lets reverse...
    path.reverse()

    aln_s1 = []
    aln_s2 = []
    s1_idx = 0
    s2_idx = 0

    for p in path:
        if p == 1:
            aln_s1.append(s1[s1_idx])
            aln_s2.append(s2[s2_idx])
            s1_idx += 1
            s2_idx += 1
        elif p == 2:
            aln_s1.append("-")
            aln_s2.append(s2[s2_idx])
            s2_idx += 1
        else:
            aln_s1.append(s1[s1_idx])
            aln_s2.append("-")
            s1_idx += 1

    aln_s1 = "".join(aln_s1)
    aln_s2 = "".join(aln_s2)

    return (aln_s1, aln_s2)
