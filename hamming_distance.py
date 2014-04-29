#!/usr/bin/env python
"""Calculate the hamming distance for a list of indexes"""

__author__ = "Aaron Berlin"
__copyright__ = "Copyright 2014, Enzymatics"
__credits__ = ["Aaron Berlin"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Aaron Berlin"
__email__ = "aberlin@enzymatics.com"
__status__ = "Development"


from itertools import imap
import operator
import sys


def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2))


def main(args):

    index_list = []
    with open(args[1]) as infile:
        for line in infile.readlines():
            index_list.append(line.rstrip())
        for i, index in enumerate(index_list):
            min_distance = 10
            total = 0
            count = 0
            for j, comp in enumerate(index_list):
                if not i == j:
                    distance = hamming_distance(index, comp)
                    if distance < min_distance:
                        min_distance = distance
                    count += 1
                    total += distance
            print index, total / count, min_distance


if __name__ == "__main__":
    main(sys.argv)