#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == "PYQT5":
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == "PYQT4":
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception("Unsupported Version of PyQt: {}".format(PYQT_VER))

import math
import time

from unrestricted import unrestricted_alignment
from banded import banded_alignment


class GeneSequencing:

    def __init__(self):
        pass

    # This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
    # handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    def align(self, sequences, table, banded, align_length):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):

                if j < i:
                    s = {}
                else:
                    ###################################################################################################
                    # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
                    # score = i+j;
                    # alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                    # 	len(sequences[i]), align_length, ',BANDED' if banded else '')
                    # alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                    # 	len(sequences[j]), align_length, ',BANDED' if banded else '')
                    if banded:
                        score, alignment1, alignment2 = banded_alignment(
                            sequences[i], sequences[j], align_length, banded
                        )
                    else:
                        score, alignment1, alignment2 = unrestricted_alignment(
                            sequences[i], sequences[j], align_length
                        )
                    ###################################################################################################
                    s = {
                        "align_cost": score,
                        "seqi_first100": alignment1,
                        "seqj_first100": alignment2,
                    }

                    table.item(i, j).setText(
                        "{}".format(int(score) if score != math.inf else score)
                    )
                    table.update()
                jresults.append(s)
            results.append(jresults)
        return results
