from typing import Tuple

def Backtrack(backtrack, scores, v, w, i, j, gap):
    seq1 = []
    seq2 = []
    result = scores[len(v)][len(w)]
    while i > 0 or j > 0:
        if backtrack[i][j] == "diagonal":
            seq1.append(v[i - 1])
            seq2.append(w[j - 1])
            i -= 1
            j -= 1
        elif j == 0 or backtrack[i][j] == "down":
            seq1.append(v[i - 1])
            seq2.append(gap)
            i -= 1
        else:
            seq1.append(gap)
            seq2.append(w[j - 1])
            j -= 1
    seq1 = "".join(seq1[::-1])
    seq2 = "".join(seq2[::-1])
    return result, seq1, seq2

def unrestricted_alignment(match_reward: int, mismatch_penalty: int, indel_penalty: int,
                     s: str, t: str, gap) -> Tuple[int, str, str]:
    backtrack = []
    table = []
    for u in range(len(s)+1):
        row = []
        back = []
        for y in range(len(t)+1):
            if u == 0:
                row.append(indel_penalty*y)
                back.append('')
            elif y == 0:
                row.append(indel_penalty*u)
                back.append('')
            else:
                row.append(0)
                back.append('')
        table.append(row)
        backtrack.append(back)
    i = 1
    while i < len(s)+1:
        j = 1
        while j < len(t)+1:
            match = mismatch_penalty
            if s[i-1] == t[j-1]:
                match = match_reward
            table[i][j] = min(table[i-1][j-1]+match, table[i-1][j]+indel_penalty, table[i][j-1]+indel_penalty)
            if table[i][j] == table[i-1][j-1]+match:
                backtrack[i][j] = "diagonal"
            elif table[i][j] == table[i][j-1]+indel_penalty:
                backtrack[i][j] = "right"
            elif table[i][j] == table[i-1][j]+indel_penalty:
                backtrack[i][j] = "down"
            j+=1
        i+=1
    result1 = Backtrack(backtrack, table, s, t, len(s), len(t), gap)
    return result1

def BacktrackBanded(backtrack, scores, v, w, i, j, gap, starts):
    seq1 = []
    seq2 = []
    result = scores[i][j]
    while i > 0 or j > 0:
        starting_point = starts[i]
        adjusted_j = j + starting_point
        if backtrack[i][j] == "diagonal":
            seq1.append(v[i - 1])
            seq2.append(w[adjusted_j - 1])
            i -= 1
            j = adjusted_j - starts[i] - 1
        elif backtrack[i][j] == "right":
            seq1.append(gap)
            seq2.append(w[adjusted_j - 1])
            j -= 1
        elif backtrack[i][j] == "down":
            seq1.append(v[i - 1])
            seq2.append(gap)
            i -= 1
            j = adjusted_j - starts[i]
    seq1 = "".join(seq1[::-1])
    seq2 = "".join(seq2[::-1])
    return result, seq1, seq2

def banded_aligntment(match_award, sub_penalty, indel_penalty, banded_width, seq1, seq2, gap):
    width = banded_width * 2 + 1
    backtrack = []
    table = []
    starts = []
    for u in range(len(seq1) + 1):
        row = []
        back = []
        begin = u - banded_width
        row_start = max(0, begin)
        row_end = len(seq2) + 1
        if begin <= 0:
            row_end = min(len(seq2) + 1, width - abs(begin))
        else:
            row_end = min(len(seq2) + 1, begin + width)
        starts.append(row_start)
        for y in range(row_start, row_end):
            if u == 0:
                row.append(indel_penalty * y)
                back.append('right')
            elif y == 0:
                row.append(indel_penalty * u)
                back.append('down')
            else:
                row.append(0)
                back.append('')
        table.append(row)
        backtrack.append(back)
    i = 1
    while i < len(seq1) + 1:
        starting_point = starts[i]
        if starting_point == 0:
            temp = 1
        else:
            temp = 0
        while temp < len(table[i]):
            j = temp + starting_point
            match = sub_penalty
            if seq1[i-1] == seq2[j-1]:
                match = match_award
            j = min(j, starts[i] + len(table[i]) - 1)
            diag = float('inf')
            up = float('inf')
            left = float('inf')
            if starts[i-1] <= j - 1 < starts[i-1] + len(table[i-1]):
                diag = table[i-1][j-1 - starts[i-1]] + match
            if starts[i-1] <= j < starts[i-1] + len(table[i-1]):
                up = table[i-1][j - starts[i-1]] + indel_penalty
            if starts[i] <= j - 1 < starts[i] + len(table[i]):
                left = table[i][temp - 1] + indel_penalty
            table[i][temp] = min(diag, up, left)
            if table[i][temp] == diag:
                backtrack[i][temp] = "diagonal"
            elif table[i][temp] == left:
                backtrack[i][temp] = "right"
            elif table[i][temp] == up:
                backtrack[i][temp] = "down"
            temp += 1
        i += 1
    v = len(table[-1]) - 1
    result1 = BacktrackBanded(backtrack, table, seq1, seq2, len(seq1), v, gap, starts)
    return result1

def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap='-'
) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap: the character to use to represent gaps in the alignment strings
        :return: alignment cost, alignment 1, alignment 2
    """
    if banded_width == -1:
        return unrestricted_alignment(match_award, sub_penalty, indel_penalty,  seq1, seq2, gap)
    else:
        return banded_aligntment(match_award, sub_penalty, indel_penalty, banded_width, seq1, seq2, gap)

match = -3
mismatch = 1 
indel = 5
s = "ATCTGG"
t = "ATGGG"

# Score: -6
# ATCTGG
# AT-GGG

print(align(s, t, match, indel, mismatch, 3))