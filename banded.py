def banded_alignment(seq1, seq2, align_length, banded):
    # Used to compute the bandwidth for banded version
    MAXINDELS = 3

    # Used to implement Needleman-Wunsch scoring
    MATCH = -3
    INDEL = 5
    SUB = 1

    seq1 = seq1[:align_length]
    seq2 = seq2[:align_length]

    seq1_len = len(seq1)
    seq2_len = len(seq2)

    def check_boundaries(row, col):
        if row < 0 or row > seq1_len or col < 0 or col > seq2_len:
            return False
        if banded == True:
            # check if row and col is within the banded region
            if col < row - MAXINDELS or col > row + MAXINDELS:
                return False
            if row < col - MAXINDELS or row > col + MAXINDELS:
                return False
        return True

    # adjusting the column
    def adjust_column(row, col):
        if row > MAXINDELS:
            return col - (row - MAXINDELS)
        return col

    # check sequence lengths and return inf if true
    if abs(seq1_len - seq2_len) > MAXINDELS:
        return float("inf"), "No Alignmnet Possible", "No Alignment Possible"

    cost = []
    # initializing cost matrix
    for row in range(seq1_len + 1):
        row_values = []
        for col in range(row - MAXINDELS, row + MAXINDELS + 1):
            if check_boundaries(row, col):
                if col == 0:
                    value = row * INDEL
                elif row == 0:
                    value = col * INDEL
                else:
                    value = float("inf")
                row_values.append(value)
        cost.append(row_values)

    backpointer = {}
    # intializing back pointer
    for row in range(seq1_len + 1):
        for col in range(seq2_len + 1):
            if col == 0:
                backpointer[(row, col)] = (row - 1, col)
            elif row == 0:
                backpointer[(row, col)] = (row, col - 1)

    # populating the cost matrix and backpointer
    for row in range(seq1_len + 1):
        for col in range(row - MAXINDELS, row + MAXINDELS + 1):
            if not check_boundaries(row, col):
                continue

            (row_curr, col_curr) = (row, adjust_column(row, col))

            if check_boundaries(row - 1, col - 1):
                cost[row_curr][col_curr] = cost[row - 1][
                    adjust_column(row - 1, col - 1)
                ] + (MATCH if seq1[row - 1] == seq2[col - 1] else SUB)
                backpointer[(row, col)] = (row - 1, col - 1)

            if check_boundaries(row - 1, col):
                prev_above_cost = cost[row - 1][adjust_column(row - 1, col)]
                new_above_cost = prev_above_cost + INDEL
                if new_above_cost <= cost[row_curr][col_curr]:
                    cost[row_curr][col_curr] = new_above_cost
                    backpointer[(row, col)] = (row - 1, col)

            if check_boundaries(row, col - 1):
                prev_left_cost = cost[row][adjust_column(row, col - 1)]
                new_left_cost = prev_left_cost + INDEL
                if new_left_cost <= cost[row_curr][col_curr]:
                    cost[row_curr][col_curr] = new_left_cost
                    backpointer[(row, col)] = (row, col - 1)

    alignment1 = ""
    alignment2 = ""
    (cur_row, cur_col) = (seq1_len, seq2_len)
    while cur_row > 0 or cur_col > 0:
        (prev_row, prev_col) = backpointer[(cur_row, cur_col)]
        if prev_row < cur_row and prev_col < cur_col:
            alignment1 += seq1[cur_row - 1]
            alignment2 += seq2[cur_col - 1]
        elif prev_row < cur_row:
            alignment1 += seq1[cur_row - 1]
            alignment2 += "-"
        elif prev_col < cur_col:
            alignment1 += "-"
            alignment2 += seq2[cur_col - 1]
        (cur_row, cur_col) = (prev_row, prev_col)

    alignment1 = alignment1[::-1][:100]
    alignment2 = alignment2[::-1][:100]

    score = cost[-1][-1]
    return score, alignment1, alignment2
