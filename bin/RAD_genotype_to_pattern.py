#!/usr/bin/env python
'''
Created on Mar 9, 2011

@author: tcezard
'''
import re
import sys
import logging
import itertools
import operator
from optparse import OptionParser
import os
from collections import Counter

from utils import utils_logging
from IO_interface import vcfIO


OUTPUT_TYPE_JOINMAP = "joinmap"
OUTPUT_TYPE_CARTHAGENE = "carthagene"
OUTPUT_TYPE_ONEMAP = "onemap"
mendelian_error = 50
OUTPUT_TYPE = [OUTPUT_TYPE_JOINMAP, OUTPUT_TYPE_CARTHAGENE, OUTPUT_TYPE_ONEMAP]

all_possible_patterns = {"0/0 0/0": {"./0": 0, "0/0": 0, "./.": 0},
                         "0/0 0/1": {"./0": 0, "0/1": 0, "0/0": 0, "./.": 0, "./1": 0},
                         "0/0 0/2": {"0/2": 0, "./0": 0, "0/0": 0, "./2": 0, "./.": 0},
                         "0/0 0/3": {"./0": 0, "./3": 0, "0/0": 0, "./.": 0, "0/3": 0},
                         "0/0 1/1": {"./0": 0, "0/1": 0, "./.": 0, "./1": 0},
                         "0/0 1/2": {"./0": 0, "./.": 0, "0/2": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "0/0 1/3": {"0/3": 0, "./.": 0, "./0": 0, "./1": 0, "0/1": 0, "./3": 0},
                         "0/0 2/2": {"0/2": 0, "./0": 0, "./2": 0, "./.": 0},
                         "0/0 2/3": {"./0": 0, "./.": 0, "0/2": 0, "0/3": 0, "./3": 0, "./2": 0},
                         "0/0 3/3": {"./0": 0, "./3": 0, "./.": 0, "0/3": 0},
                         "0/1 0/0": {"./0": 0, "0/1": 0, "0/0": 0, "./.": 0, "./1": 0},
                         "0/1 0/1": {"./.": 0, "./0": 0, "1/1": 0, "./1": 0, "0/1": 0, "0/0": 0},
                         "0/1 0/2": {"./0": 0, "./.": 0, "1/2": 0, "0/2": 0, "./1": 0, "0/1": 0, "0/0": 0, "./2": 0},
                         "0/1 0/3": {"1/3": 0, "./1": 0, "./.": 0, "./0": 0, "0/3": 0, "0/1": 0, "0/0": 0, "./3": 0},
                         "0/1 1/1": {"./0": 0, "0/1": 0, "1/1": 0, "./.": 0, "./1": 0},
                         "0/1 1/2": {"./0": 0, "./.": 0, "1/2": 0, "0/2": 0, "1/1": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "0/1 1/3": {"1/3": 0, "./1": 0, "./.": 0, "./0": 0, "1/1": 0, "0/3": 0, "0/1": 0, "./3": 0},
                         "0/1 2/2": {"./0": 0, "./.": 0, "1/2": 0, "0/2": 0, "./1": 0, "./2": 0},
                         "0/1 2/3": {"./0": 0, "1/3": 0, "./1": 0, "./.": 0, "1/2": 0, "0/2": 0, "0/3": 0, "./3": 0,
                                     "./2": 0},
                         "0/1 3/3": {"1/3": 0, "./1": 0, "./.": 0, "./0": 0, "0/3": 0, "./3": 0},
                         "0/2 0/0": {"./0": 0, "0/2": 0, "0/0": 0, "./2": 0, "./.": 0},
                         "0/2 0/1": {"0/2": 0, "./.": 0, "1/2": 0, "./0": 0, "./1": 0, "0/1": 0, "0/0": 0, "./2": 0},
                         "0/2 0/2": {"./0": 0, "2/2": 0, "0/2": 0, "./.": 0, "0/0": 0, "./2": 0},
                         "0/2 0/3": {"0/2": 0, "./.": 0, "2/3": 0, "./0": 0, "0/3": 0, "./3": 0, "0/0": 0, "./2": 0},
                         "0/2 1/1": {"./.": 0, "1/2": 0, "./0": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "0/2 1/2": {"./0": 0, "2/2": 0, "1/2": 0, "0/2": 0, "./.": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "0/2 1/3": {"./1": 0, "./.": 0, "1/2": 0, "./0": 0, "0/3": 0, "0/1": 0, "./2": 0, "2/3": 0,
                                     "./3": 0},
                         "0/2 2/2": {"0/2": 0, "./0": 0, "./2": 0, "./.": 0, "2/2": 0},
                         "0/2 2/3": {"./0": 0, "2/2": 0, "0/2": 0, "./.": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "0/2 3/3": {"./.": 0, "./0": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "0/3 0/0": {"./0": 0, "./3": 0, "0/0": 0, "./.": 0, "0/3": 0},
                         "0/3 0/1": {"1/3": 0, "0/3": 0, "./1": 0, "./.": 0, "./0": 0, "./3": 0, "0/1": 0, "0/0": 0},
                         "0/3 0/2": {"./2": 0, "./0": 0, "./.": 0, "0/2": 0, "0/3": 0, "./3": 0, "0/0": 0, "2/3": 0},
                         "0/3 0/3": {"3/3": 0, "./.": 0, "./0": 0, "0/3": 0, "./3": 0, "0/0": 0},
                         "0/3 1/1": {"1/3": 0, "./.": 0, "./0": 0, "./1": 0, "0/1": 0, "./3": 0},
                         "0/3 1/2": {"./0": 0, "1/3": 0, "./1": 0, "./.": 0, "0/2": 0, "./3": 0, "0/1": 0, "./2": 0,
                                     "2/3": 0},
                         "0/3 1/3": {"3/3": 0, "1/3": 0, "0/3": 0, "./1": 0, "./.": 0, "./0": 0, "./3": 0, "0/1": 0},
                         "0/3 2/2": {"./0": 0, "./.": 0, "0/2": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "0/3 2/3": {"./0": 0, "3/3": 0, "./.": 0, "0/2": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "0/3 3/3": {"./0": 0, "./3": 0, "./.": 0, "3/3": 0, "0/3": 0},
                         "1/1 0/0": {"./0": 0, "0/1": 0, "./.": 0, "./1": 0},
                         "1/1 0/1": {"./0": 0, "0/1": 0, "1/1": 0, "./.": 0, "./1": 0},
                         "1/1 0/2": {"./.": 0, "1/2": 0, "./0": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "1/1 0/3": {"1/3": 0, "./1": 0, "./.": 0, "./0": 0, "./3": 0, "0/1": 0},
                         "1/1 1/1": {"1/1": 0, "./.": 0, "./1": 0},
                         "1/1 1/2": {"1/1": 0, "./1": 0, "./2": 0, "./.": 0, "1/2": 0},
                         "1/1 1/3": {"1/3": 0, "./3": 0, "1/1": 0, "./.": 0, "./1": 0},
                         "1/1 2/2": {"./2": 0, "./1": 0, "./.": 0, "1/2": 0},
                         "1/1 2/3": {"1/3": 0, "./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0},
                         "1/1 3/3": {"1/3": 0, "./3": 0, "./.": 0, "./1": 0},
                         "1/2 0/0": {"./0": 0, "./.": 0, "0/2": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "1/2 0/1": {"./0": 0, "./.": 0, "1/2": 0, "0/2": 0, "1/1": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "1/2 0/2": {"./0": 0, "2/2": 0, "1/2": 0, "0/2": 0, "./.": 0, "./1": 0, "0/1": 0, "./2": 0},
                         "1/2 0/3": {"./0": 0, "1/3": 0, "./.": 0, "0/2": 0, "./1": 0, "0/1": 0, "./2": 0, "2/3": 0,
                                     "./3": 0},
                         "1/2 1/1": {"1/1": 0, "1/2": 0, "./2": 0, "./.": 0, "./1": 0},
                         "1/2 1/2": {"./.": 0, "1/2": 0, "1/1": 0, "2/2": 0, "./1": 0, "./2": 0},
                         "1/2 1/3": {"1/3": 0, "./.": 0, "1/2": 0, "1/1": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/2 2/2": {"./2": 0, "./1": 0, "./.": 0, "2/2": 0, "1/2": 0},
                         "1/2 2/3": {"1/3": 0, "2/2": 0, "1/2": 0, "./.": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/2 3/3": {"1/3": 0, "./.": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/3 0/0": {"./1": 0, "./.": 0, "./0": 0, "0/3": 0, "0/1": 0, "./3": 0},
                         "1/3 0/1": {"1/3": 0, "./1": 0, "./.": 0, "./0": 0, "1/1": 0, "0/3": 0, "0/1": 0, "./3": 0},
                         "1/3 0/2": {"./1": 0, "./.": 0, "1/2": 0, "./0": 0, "0/3": 0, "0/1": 0, "./2": 0, "2/3": 0,
                                     "./3": 0},
                         "1/3 0/3": {"3/3": 0, "1/3": 0, "./1": 0, "./.": 0, "./0": 0, "0/3": 0, "0/1": 0, "./3": 0},
                         "1/3 1/1": {"1/3": 0, "./3": 0, "1/1": 0, "./.": 0, "./1": 0},
                         "1/3 1/2": {"1/3": 0, "./.": 0, "1/2": 0, "1/1": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/3 1/3": {"3/3": 0, "1/3": 0, "./.": 0, "1/1": 0, "./1": 0, "./3": 0},
                         "1/3 2/2": {"./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/3 2/3": {"3/3": 0, "1/3": 0, "./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "1/3 3/3": {"1/3": 0, "./3": 0, "./.": 0, "3/3": 0, "./1": 0},
                         "2/2 0/0": {"0/2": 0, "./0": 0, "./2": 0, "./.": 0},
                         "2/2 0/1": {"./0": 0, "./.": 0, "1/2": 0, "0/2": 0, "./1": 0, "./2": 0},
                         "2/2 0/2": {"0/2": 0, "./0": 0, "./2": 0, "./.": 0, "2/2": 0},
                         "2/2 0/3": {"./0": 0, "./.": 0, "0/2": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/2 1/1": {"./2": 0, "./1": 0, "./.": 0, "1/2": 0},
                         "2/2 1/2": {"./2": 0, "./1": 0, "./.": 0, "2/2": 0, "1/2": 0},
                         "2/2 1/3": {"./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/2 2/2": {"./2": 0, "./.": 0, "2/2": 0},
                         "2/2 2/3": {"./.": 0, "./3": 0, "./2": 0, "2/3": 0, "2/2": 0},
                         "2/2 3/3": {"./3": 0, "./2": 0, "2/3": 0, "./.": 0},
                         "2/3 0/0": {"./0": 0, "./.": 0, "0/2": 0, "0/3": 0, "./3": 0, "./2": 0},
                         "2/3 0/1": {"./0": 0, "1/3": 0, "./1": 0, "./.": 0, "1/2": 0, "0/2": 0, "0/3": 0, "./3": 0,
                                     "./2": 0},
                         "2/3 0/2": {"./0": 0, "2/2": 0, "0/2": 0, "./.": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/3 0/3": {"./0": 0, "3/3": 0, "./.": 0, "0/2": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/3 1/1": {"1/3": 0, "./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0},
                         "2/3 1/2": {"1/3": 0, "2/2": 0, "1/2": 0, "./.": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/3 1/3": {"3/3": 0, "1/3": 0, "./.": 0, "1/2": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/3 2/2": {"./.": 0, "./3": 0, "./2": 0, "2/3": 0, "2/2": 0},
                         "2/3 2/3": {"3/3": 0, "2/2": 0, "./.": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "2/3 3/3": {"./.": 0, "./3": 0, "./2": 0, "2/3": 0, "3/3": 0},
                         "3/3 0/0": {"./0": 0, "./3": 0, "./.": 0, "0/3": 0},
                         "3/3 0/1": {"1/3": 0, "0/3": 0, "./.": 0, "./0": 0, "./1": 0, "./3": 0},
                         "3/3 0/2": {"./.": 0, "./0": 0, "0/3": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "3/3 0/3": {"./0": 0, "./3": 0, "./.": 0, "3/3": 0, "0/3": 0},
                         "3/3 1/1": {"1/3": 0, "./3": 0, "./.": 0, "./1": 0},
                         "3/3 1/2": {"1/3": 0, "./.": 0, "./1": 0, "./3": 0, "./2": 0, "2/3": 0},
                         "3/3 1/3": {"1/3": 0, "./3": 0, "./.": 0, "3/3": 0, "./1": 0},
                         "3/3 2/2": {"./3": 0, "./2": 0, "2/3": 0, "./.": 0},
                         "3/3 2/3": {"./.": 0, "./3": 0, "./2": 0, "2/3": 0, "3/3": 0},
                         "3/3 3/3": {"./3": 0, "./.": 0, "3/3": 0}}


def generate_empty_hash_with_sample(all_samples):
    hash = {}
    for sample in all_samples:
        hash[sample] = {}
    hash['all'] = {}
    return hash


def count_with_hash(hashtable, key):
    if hashtable.has_key(key):
        hashtable[key] += 1
    else:
        hashtable[key] = 1


def get_offsprings_possible_genotypes_jm():
    # the parents pattern assume that the mother is first and the father second
    jm_parents_pattern = {"0/1 0/0": {"type": "<lmxll>", "0/0": "ll", "0/1": "lm", "1/0": "lm"},
                          "0/1 1/1": {"type": "<lmxll>", "1/1": "ll", "0/1": "lm", "1/0": "lm"},
                          "0/0 0/1": {"type": "<nnxnp>", "0/0": "nn", "0/1": "np", "1/0": "np"},
                          "1/1 0/1": {"type": "<nnxnp>", "1/1": "nn", "0/1": "np", "1/0": "np"},
                          "0/1 0/1": {"type": "<hkxhk>", "0/0": "hh", "0/1": "hk", "1/0": "hk", "1/1": "kk"},
                          "0/1 0/2": {"type": "<efxeg>", "0/0": "ee", "0/1": "ef", "0/2": "eg", "1/2": "fg"},
                          "0/1 1/2": {"type": "<efxeg>", "1/1": "ee", "0/1": "ef", "0/2": "fg", "1/2": "eg"},
                          "0/2 0/1": {"type": "<efxeg>", "0/0": "ee", "0/1": "eg", "0/2": "ef", "1/2": "fg"},
                          "0/2 1/2": {"type": "<efxeg>", "2/2": "ee", "0/1": "fg", "0/2": "ef", "1/2": "eg"},
                          "0/1 2/3": {"type": "<abxcd>", "0/2": "ac", "0/3": "ad", "1/2": "bc", "1/3": "bd"},
                          "0/2 1/3": {"type": "<abxcd>", "0/1": "ac", "0/3": "ad", "1/2": "bc", "2/3": "bd"},
                          "0/3 1/2": {"type": "<abxcd>", "0/1": "ac", "0/2": "ad", "1/3": "bc", "2/3": "bd"},
                          "1/2 0/3": {"type": "<abxcd>", "0/1": "ac", "1/3": "ad", "0/2": "bc", "2/3": "bd"},
                          "1/3 0/2": {"type": "<abxcd>", "0/1": "ac", "0/2": "ad", "0/3": "bc", "2/3": "bd"},
                          "2/3 0/1": {"type": "<abxcd>", "0/2": "ac", "1/2": "ad", "0/3": "bc", "1/3": "bd"},
                          #"0/0 1/1":{},
                          #"0/0 1/2":{},
                          #"1/1 0/0":{},
                          #"1/2 0/0":{},
    }

    return jm_parents_pattern


def get_offsprings_possible_genotypes_om():
    # the parents pattern assume that the mother is first and the father second

    om_parents_pattern = {"0/1 0/0": {"type": "D1.10", "0/0": "a", "0/1": "ab", "1/0": "ab"},
                          "0/1 1/1": {"type": "D1.10", "1/1": "a", "0/1": "ab", "1/0": "ab"},
                          "0/0 0/1": {"type": "D1.10", "0/0": "a", "0/1": "ab", "1/0": "ab"},
                          "1/1 0/1": {"type": "D1.10", "1/1": "a", "0/1": "ab", "1/0": "ab"},
                          "0/1 0/1": {"type": "B3.7", "0/0": "a", "0/1": "ab", "1/0": "ab", "1/1": "b"},
                          "0/1 0/2": {"type": "A.2", "0/0": "a", "0/1": "ba", "0/2": "ac", "1/2": "bc"},
                          "0/1 1/2": {"type": "A.2", "1/1": "a", "0/1": "ba", "0/2": "bc", "1/2": "ac"},
                          "0/2 0/1": {"type": "A.2", "0/0": "a", "0/1": "ac", "0/2": "ba", "1/2": "bc"},
                          "0/2 1/2": {"type": "A.2", "2/2": "a", "0/1": "bc", "0/2": "ba", "1/2": "ac"},
                          "0/1 2/3": {"type": "A.1", "0/2": "ac", "0/3": "ad", "1/2": "bc", "1/3": "bd"},
                          "0/2 1/3": {"type": "A.1", "0/1": "ac", "0/3": "ad", "1/2": "bc", "2/3": "bd"},
                          "0/3 1/2": {"type": "A.1", "0/1": "ac", "0/2": "ad", "1/3": "bc", "2/3": "bd"},
                          "1/2 0/3": {"type": "A.1", "0/1": "ac", "1/3": "ad", "0/2": "bc", "2/3": "bd"},
                          "1/3 0/2": {"type": "A.1", "0/1": "ac", "0/2": "ad", "0/3": "bc", "2/3": "bd"},
                          "2/3 0/1": {"type": "A.1", "0/2": "ac", "1/2": "ad", "0/3": "bc", "1/3": "bd"},
                          "0/0 1/2": {"type": "D2.14", "0/1": "ac", "0/2": "bc"},
                          "1/2 0/0": {"type": "D2.14", "0/1": "ac", "0/2": "bc"}
    }

    return om_parents_pattern


def get_offsprings_possible_genotypes_ca():
    # the parents pattern assume that the mother is first and the father second

    female = {
        #"0/0 0/1":{},                          Not informative in the female
        #"1/1 0/1":{},                          Not informative in the female
        "0/1 0/0": {"0/0": "A", "0/1": "H"},
        "0/1 1/1": {"1/1": "H", "0/1": "A"},

        "0/1 0/1": {"0/0": "A", "0/1": "-", "1/1": "H"},

        "0/1 0/2": {"0/0": "A", "0/1": "H", "0/2": "A", "1/2": "H"},
        "0/1 1/2": {"0/1": "A", "0/2": "A", "1/1": "H", "1/2": "H"},
        "0/2 0/1": {"0/0": "A", "0/1": "A", "0/2": "H", "1/2": "H"},
        "0/2 1/2": {"0/1": "A", "0/2": "A", "1/2": "H", "1/2": "H"},
        "1/2 0/1": {"0/1": "A", "1/1": "A", "0/2": "H", "1/2": "H"},
        "1/2 0/2": {"0/1": "A", "1/2": "A", "0/2": "H", "2/2": "H"},

        "0/1 2/3": {"0/2": "A", "0/3": "A", "1/2": "H", "1/3": "H"},
        "0/2 1/3": {"0/1": "A", "0/3": "A", "1/2": "H", "2/3": "H"},
        "0/3 1/2": {"0/1": "A", "0/2": "A", "1/3": "H", "2/3": "H"},
        "1/2 0/3": {"0/1": "A", "1/3": "A", "0/2": "H", "2/3": "H"},
        "1/3 0/2": {"0/1": "A", "1/2": "A", "0/3": "H", "2/3": "H"},
        "2/3 0/1": {"0/2": "A", "1/2": "A", "0/3": "H", "1/3": "H"},

        "1/2 0/0": {"0/1": "A", "0/1": "H"}
        #"0/0 1/2":{},                           Not informative in the female
    }

    male = {
        "0/0 0/1": {"0/0": "A", "0/1": "H"},
        "1/1 0/1": {"0/1": "A", "1/1": "H"},
        #"0/1 0/0":                              Not informative in the male
        #"0/1 1/1":                              Not informative in the male

        "0/1 0/1": {"0/0": "A", "0/1": "-", "1/1": "H"},

        "0/1 0/2": {"0/0": "A", "0/1": "A", "0/2": "H", "1/2": "H"},
        "0/1 1/2": {"0/1": "A", "1/1": "A", "0/2": "H", "1/2": "H"},
        "0/2 0/1": {"0/0": "A", "0/2": "A", "0/1": "H", "1/2": "H"},
        "0/2 1/2": {"0/1": "A", "1/2": "A", "0/2": "H", "1/2": "H"},
        "1/2 0/1": {"0/1": "A", "0/2": "A", "1/1": "H", "1/2": "H"},
        "1/2 0/2": {"0/1": "A", "0/2": "A", "1/2": "H", "2/2": "H"},

        "0/1 2/3": {"0/2": "A", "1/2": "A", "0/3": "H", "1/3": "H"},
        "0/2 1/3": {"0/1": "A", "1/2": "A", "0/3": "H", "2/3": "H"},
        "0/3 1/2": {"0/1": "A", "1/3": "A", "0/2": "H", "2/3": "H"},
        "1/2 0/3": {"0/1": "A", "0/2": "A", "1/3": "H", "2/3": "H"},
        "1/3 0/2": {"0/1": "A", "0/3": "A", "1/2": "H", "2/3": "H"},
        "2/3 0/1": {"0/2": "A", "0/3": "A", "1/2": "H", "1/3": "H"},

        #"1/2 0/0":{}                              Not informative in the male
        "0/0 1/2": {"0/1": "A", "0/1": "H"}
    }
    return female, male


def load_sex_info(file):
    sample_to_sex = {}
    sex_to_sample = {}
    sex_to_sample["M"] = []
    sex_to_sample["F"] = []
    ordered_samples = []
    open_file = utils_logging.open_input_file(file, pipe=False)
    for line in open_file:
        sp_line = line.strip().split()
        ordered_samples.append(sp_line[0])
        if sp_line[1] == "m" or sp_line[1] == "M":
            sex_to_sample["M"].append(sp_line[0])
            sample_to_sex[sp_line[0]] = "M"
        elif sp_line[1] == "f" or sp_line[1] == "F":
            sex_to_sample["F"].append(sp_line[0])
            sample_to_sex[sp_line[0]] = "F"

    open_file.close()
    return sample_to_sex, sex_to_sample, ordered_samples


def reverse_consensus(consensus):
    tmp = []
    for letter in consensus:
        if letter == 'A':
            tmp.append('H')
        elif letter == 'H':
            tmp.append('A')
        else:
            tmp.append(letter)
    return "".join(tmp)


class Phased_vcf_reader():
    def __init__(self, vcf_file, mother, father, geno_qual_threshold):

        self.file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
        self.reader = vcfIO.VcfReader(self.file_handle)
        self.geno_qual_threshold = geno_qual_threshold
        self.mother = mother
        self.father = father
        self.nb_snps = 0
        self.marker_created = 0
        self.discard_quality_parent = 0
        self.ordered_sample = []
        for sample in self.reader.get_sample_names():
            self.ordered_sample.append(sample)
        self.ordered_sample.remove(self.mother)
        self.ordered_sample.remove(self.father)

    def read(self):

        end_of_phase = False
        phased_genotype_mother = []
        phased_genotype_father = []
        phased_vcf_record = []
        previous_reference = None
        start_session_mother = False
        start_session_father = False

        for vcf_records in self.reader:
            # First check that the parent are callable
            self.nb_snps += 1
            discarded_snp = False
            geno_mother = vcf_records.get_genotype(self.mother)
            gq_mother = vcf_records.get_genotype_quality(self.mother)
            geno_father = vcf_records.get_genotype(self.father)
            gq_father = vcf_records.get_genotype_quality(self.father)
            if gq_mother < self.geno_qual_threshold or gq_father < self.geno_qual_threshold:
                if gq_mother < self.geno_qual_threshold and gq_father < self.geno_qual_threshold:
                    logging.info(
                        "Bad quality in both parent: %s:%s" % (vcf_records.get_reference(), vcf_records.get_position()))
                elif gq_father < self.geno_qual_threshold:
                    logging.info(
                        "Bad quality in Father: %s:%s" % (vcf_records.get_reference(), vcf_records.get_position()))
                elif gq_mother < self.geno_qual_threshold:
                    logging.info(
                        "Bad quality in Mother: %s:%s" % (vcf_records.get_reference(), vcf_records.get_position()))
                self.discard_quality_parent += 1
                discarded_snp = True
                self.marker_created += 1
                #break the Phase if we the quality in the parent is bad
                #do not add this SNP
                add_SNP_before = False
                add_SNP_after = False
                end_of_phase = True
                #logging.debug("stop phasing because %s:%s is removed"%(vcf_records.get_reference(),vcf_records.get_position()))
            elif vcf_records.is_phased_with_previous(sample=self.mother) and vcf_records.is_phased_with_previous(
                    sample=self.father):
                #Keep the Phase if both parents are keeping the phase
                #Add this SNP to the phased set
                add_SNP_before = True
                add_SNP_after = False
                end_of_phase = False

            else:
                #Break the Phase if at least one parent is breaking the phase
                #Add this SNP to a new phased set
                add_SNP_before = False
                add_SNP_after = True
                #to correctly phase homozygous location at the beginning of contigs 
                if not vcf_records.is_phased_with_previous(sample=self.mother):
                    start_session_mother = True
                if not vcf_records.is_phased_with_previous(sample=self.father):
                    start_session_father = True
                if start_session_mother and start_session_father:
                    #logging.debug("stop phasing because not phased with previous %s:%s"%(vcf_records.get_reference(), vcf_records.get_position()))
                    end_of_phase = True
                else:
                    end_of_phase = False

            if previous_reference and previous_reference != vcf_records.get_reference():
                #break the Phase if we change reference
                add_SNP_before = False
                if not discarded_snp:
                    add_SNP_after = True
                end_of_phase = True
                #logging.debug("stop phasing because changing reference: now %s --> %s"%(previous_reference, vcf_records.get_reference()))
                start_session_mother = start_session_father = False

            if vcf_records.get_info_tag_value('PhasingInconsistent'):
                #break the phase if GATK says the phasing is inconsistent
                #logging.debug("stop phasing because Found a PhasingInconsistent mark: %s:%s"%(vcf_records.get_reference(), vcf_records.get_position()))
                end_of_phase = True

            if add_SNP_before:
                phased_genotype_mother.append(re.split(r'[|/]', geno_mother))
                phased_genotype_father.append(re.split(r'[|/]', geno_father))
                phased_vcf_record.append(vcf_records)
            #print previous_reference, vcf_records.get_reference() , vcf_records.get_position()
            #print 'in phase for mother=%s, in phase for father=%s, end_of_phase=%s, add_SNP_before=%s, add_SNP_after=%s'%(
            #                                                         vcf_records.is_phased_with_previous(sample=self.mother),
            #                                                         vcf_records.is_phased_with_previous(sample=self.father),
            #                                                         end_of_phase, add_SNP_before,add_SNP_after)
            #print 'start_session_mother=%s, start_session_father=%s'%(start_session_mother,start_session_father)
            #print 'nb_of SNPs in phased set = %s'%len(phased_vcf_record)
            #if  len(phased_vcf_record)>0:
            #    print ' || '.join(['%s--%s'%(tmp.get_reference(),tmp.get_position()) for tmp in phased_vcf_record])
            if end_of_phase and len(phased_vcf_record) > 0:
                self.marker_created += 1
                logging.debug("Create marker: %s" % (get_list_of_SNPs_name(phased_vcf_record)))
                yield (phased_genotype_mother, phased_genotype_father, phased_vcf_record)

                phased_genotype_mother = []
                phased_genotype_father = []
                phased_vcf_record = []

            if add_SNP_after:
                phased_genotype_mother.append(re.split(r'[|/]', geno_mother))
                phased_genotype_father.append(re.split(r'[|/]', geno_father))
                phased_vcf_record.append(vcf_records)
            previous_reference = vcf_records.get_reference()
        if len(phased_vcf_record) > 0:
            self.marker_created += 1
            yield (phased_genotype_mother, phased_genotype_father, phased_vcf_record)


    def get_ordered_samples(self):
        return self.ordered_sample

    def get_nb_marker_created(self):
        return self.marker_created

    def get_discard_quality_parent(self):
        return self.discard_quality_parent

    def get_nb_snps(self):
        return self.nb_snps


def phased_snps_to_markers(vcf_file, ordered_sample, mother, father, geno_qual_threshold, max_missing_offspring,
                           nb_ME_error=1):
    phased_reader = Phased_vcf_reader(vcf_file, mother, father, geno_qual_threshold)
    # non_supported_genotypes=Counter()
    #non_supported_markers=0
    broken_down_markers = 0
    broken_down_SNPs = 0
    rescued_markers = 0
    rescued_SNPs = 0
    created_markers = 0
    droped_markers = 0
    droped_SNPs = 0
    mendelian_error_count = 0
    mendelian_error_offsprings_count = Counter()
    for phased_genotype_mother, phased_genotype_father, phased_vcf_record in phased_reader.read():
        if len(phased_vcf_record) == 1:
            #If it is not a multi loci allele can allow mendelian error to go through
            #They'll be caught at the genotype matching stage
            allow_mendelian_error = True
        else:
            allow_mendelian_error = False

        geno_parents, all_offsprings_genotypes = compare_offspring_haplotype_to_parents(phased_genotype_mother,
                                                                                        phased_genotype_father,
                                                                                        phased_vcf_record,
                                                                                        ordered_sample,
                                                                                        max_missing_offspring=max_missing_offspring,
                                                                                        allow_mendelian_error=allow_mendelian_error)
        mendelian_error_offsprings = []
        if all_offsprings_genotypes:
            mendelian_error_offsprings = test_mendelian_error(geno_parents, all_offsprings_genotypes, ordered_sample)
        if all_offsprings_genotypes is None or len(mendelian_error_offsprings) > nb_ME_error:
            if len(phased_vcf_record) > 1:
                #Trying again without the phasing information
                logging.debug("breaking down marker into %s pieces" % (len(phased_vcf_record)))
                broken_down_markers += 1
                broken_down_SNPs += len(phased_vcf_record)
                marker_rescued = False
                for i in range(len(phased_vcf_record)):
                    logging.debug("marker: %s" % (get_list_of_SNPs_name([phased_vcf_record[i]])))
                    geno_parents, all_offsprings_genotypes = compare_offspring_haplotype_to_parents(
                        [phased_genotype_mother[i]],
                        [phased_genotype_father[i]],
                        [phased_vcf_record[i]], ordered_sample,
                        max_missing_offspring=max_missing_offspring,
                        allow_mendelian_error=True)
                    if all_offsprings_genotypes is not None:
                        rescued_SNPs += 1
                        marker_rescued = True
                        mendelian_error_offsprings = test_mendelian_error(geno_parents, all_offsprings_genotypes,
                                                                          ordered_sample)
                        if len(mendelian_error_offsprings) <= nb_ME_error:
                            created_markers += 1
                            yield geno_parents, all_offsprings_genotypes, [phased_vcf_record[i]]
                        else:
                            mendelian_error_count += 1
                            for offspring_error in mendelian_error_offsprings:
                                mendelian_error_offsprings_count[offspring_error] += 1
                    else:
                        droped_SNPs += 1

                if marker_rescued:
                    rescued_markers += 1
                else:
                    droped_markers += 1
            else:
                droped_markers += 1
        else:
            if len(mendelian_error_offsprings) <= nb_ME_error:
                created_markers += 1
                yield geno_parents, all_offsprings_genotypes, phased_vcf_record
            else:
                mendelian_error_count += 1
                for offspring_error in mendelian_error_offsprings:
                    mendelian_error_offsprings_count[offspring_error] += 1
    logging.info("%s SNPs discarded because of quality issue in parents" % phased_reader.get_discard_quality_parent())
    logging.info(
        '%s marker created by phasing %s SNPs' % (phased_reader.get_nb_marker_created(), phased_reader.get_nb_snps()))
    if broken_down_markers > 0:
        logging.info('%s marker broken down in %s SNPs because too many missing genotype' % (
        broken_down_markers, broken_down_SNPs))
        logging.info('%s SNPs, %s markers rescued' % (rescued_SNPs, rescued_markers))
        logging.info('%s SNPs, %s markers more than %s missing offspring ' % (
        droped_SNPs, droped_markers, max_missing_offspring))
        logging.info('%s SNPs show more than %s mendelian errors' % (mendelian_error_count, nb_ME_error))
    else:
        logging.info('%s markers more than %s missing offspring ' % (droped_markers, max_missing_offspring))
        logging.info('%s markers show more than %s mendelian errors' % (mendelian_error_count, nb_ME_error))
    logging.info('Sample with mendelian error:\n%s' % ('\n'.join(
        ['%s : %s' % (key, mendelian_error_offsprings_count.get(key)) for key in
         mendelian_error_offsprings_count.keys()])))
    logging.info('%s marker pased on after filtering' % (created_markers))


def test_mendelian_error(geno_parents, all_offsprings_genotypes, ordered_sample):
    possible_geno_offsprings = all_possible_patterns.get(geno_parents)
    samples_with_mendelian_error = []
    for i in range(len(all_offsprings_genotypes)):
        if possible_geno_offsprings.get(all_offsprings_genotypes[i]) is None:
            logging.debug('sample %s shows mendelian error with %s for parental genotype %s' % (ordered_sample[i],
                                                                                                all_offsprings_genotypes[
                                                                                                    i],
                                                                                                geno_parents))
            samples_with_mendelian_error.append(ordered_sample[i])
    return samples_with_mendelian_error


def phased_snps_to_joinmap(vcf_file, sex_info, mother, father, geno_qual_threshold, output_markers, output_vcf,
                           max_missing_offspring):
    sample_to_sex, sex_to_sample, ordered_sample = load_sex_info(sex_info)
    if mother in ordered_sample:
        ordered_sample.remove(mother)
    if father in ordered_sample:
        ordered_sample.remove(father)

    open_marker = utils_logging.open_output_file(output_markers, pipe=False)
    all_possible_genotypes = get_offsprings_possible_genotypes_jm()
    phased_markers = phased_snps_to_markers(vcf_file, ordered_sample, mother, father, geno_qual_threshold,
                                            max_missing_offspring)

    non_supported_genotypes = Counter()
    non_supported_markers = 0
    bad_pattern = 0

    for geno_parents, all_offsprings_genotypes, phased_vcf_records in phased_markers:
        name = get_list_of_SNPs_name(phased_vcf_records)
        out = output_marker_in_joinmap(name, all_possible_genotypes, ordered_sample, geno_parents,
                                       all_offsprings_genotypes, max_missing_offspring)
        if out is None:
            logging.info("Non supported genotype in parent %s for marker %s" % (geno_parents, name))
            non_supported_markers += 1
            non_supported_genotypes[geno_parents] += 1
        elif len(out) < 1:
            bad_pattern += 1
        else:
            open_marker.write(out)
    open_marker.close()
    logging.info("%s discarded because not supported by the format (%s)" % (non_supported_markers,
                                                                            ', '.join(["%s:%s" % (genotype, count) for
                                                                                       genotype, count in
                                                                                       non_supported_genotypes.items()])))
    logging.info("%s discarded because of bad pattern" % bad_pattern)


def phased_snps_to_onemap(vcf_file, sex_info, mother, father, geno_qual_threshold, output_markers, output_vcf,
                          max_missing_offspring):
    sample_to_sex, sex_to_sample, ordered_sample = load_sex_info(sex_info)
    if mother in ordered_sample:
        ordered_sample.remove(mother)
    if father in ordered_sample:
        ordered_sample.remove(father)
    possible_genotypes = get_offsprings_possible_genotypes_om()
    non_supported_genotypes = Counter()
    non_supported_markers = 0
    bad_pattern = 0
    all_markers = []

    phased_markers = phased_snps_to_markers(vcf_file, ordered_sample, mother, father, geno_qual_threshold,
                                            max_missing_offspring)
    for geno_parents, all_offsprings_genotypes, phased_vcf_record in phased_markers:
        name = get_list_of_SNPs_name(phased_vcf_record)
        out = output_marker_in_onemap(name, possible_genotypes, ordered_sample, geno_parents, all_offsprings_genotypes,
                                      max_missing_offspring)
        if out is None:
            logging.info("Non supported genotype in parent %s for marker %s" % (geno_parents, name))
            non_supported_markers += 1
            non_supported_genotypes[geno_parents] += 1
        elif len(out) < 1:
            bad_pattern += 1
        else:
            all_markers.append(out)
    open_marker = utils_logging.open_output_file(output_markers, pipe=False)
    open_marker.write("%s %s\n" % (len(ordered_sample), len(all_markers)))
    open_marker.write("".join(all_markers))
    open_marker.close()
    print "%s discarded because not supported by the format (%s)" % (non_supported_markers,
                                                                     ', '.join(["%s:%s" % (genotype, count) for
                                                                                genotype, count in
                                                                                non_supported_genotypes.items()]))
    print "%s discarded because of bad pattern" % bad_pattern


def phased_snps_to_carthagene(vcf_file, sex_info, mother, father, geno_qual_threshold, output_markers_female,
                              output_markers_male,
                              output_vcf, max_missing_offspring):
    sample_to_sex, sex_to_sample, ordered_sample = load_sex_info(sex_info)
    if mother in ordered_sample:
        ordered_sample.remove(mother)
    if father in ordered_sample:
        ordered_sample.remove(father)
    # open_output = utils_logging.open_output_file(output_vcf, pipe=False)
    possible_genotypes = get_offsprings_possible_genotypes_ca()
    sex1 = "".join([sample_to_sex.get(sample).replace("F", "A").replace("M", "H") for sample in ordered_sample])
    sex2 = "".join([sample_to_sex.get(sample).replace("M", "A").replace("F", "H") for sample in ordered_sample])

    phased_markers = phased_snps_to_markers(vcf_file, ordered_sample, mother, father, geno_qual_threshold,
                                            max_missing_offspring)

    non_supported_genotypes = Counter()
    non_supported_markers = 0
    bad_pattern = 0
    all_markers_female = []
    all_markers_male = []
    for geno_parents, all_offsprings_genotypes, phased_vcf_record in phased_markers:
        name = get_list_of_SNPs_name(phased_vcf_record)
        output_female, output_male = output_marker_in_cartagene(name, possible_genotypes, ordered_sample, geno_parents,
                                                                all_offsprings_genotypes)
        if output_female is None and output_male is None:
            logging.info("Non supported genotype in parent %s for marker %s" % (geno_parents, name))
            non_supported_markers += 1
            non_supported_genotypes[geno_parents] += 1
        else:
            if output_female:
                all_markers_female.append(output_female)
            if output_male:
                all_markers_male.append(output_male)
    all_markers_female_str = "".join(all_markers_female)
    nb_marker_female = len(all_markers_female_str.split("\n")) + 2
    open_marker_female = utils_logging.open_output_file(output_markers_female, pipe=False)
    open_marker_female.write("data type f2 backcross\n")
    open_marker_female.write("%s %s 0 0 0\n" % (len(ordered_sample), nb_marker_female))
    open_marker_female.write("*sex\t%s\n" % (sex1))
    open_marker_female.write("*sex_rev\t%s\n" % (sex2))
    open_marker_female.write(all_markers_female_str)
    open_marker_female.close()

    all_markers_male_str = "".join(all_markers_male)
    nb_marker_male = len(all_markers_male_str.split("\n")) + 2
    open_marker_male = utils_logging.open_output_file(output_markers_male, pipe=False)
    open_marker_male.write("data type f2 backcross\n")
    open_marker_male.write("%s %s 0 0 0\n" % (len(ordered_sample), nb_marker_male))
    open_marker_male.write("*sex\t%s\n" % (sex1))
    open_marker_male.write("*sex_rev\t%s\n" % (sex2))
    open_marker_male.write(all_markers_male_str)
    open_marker_male.close()

    print '%s female markers %s male markers output' % (nb_marker_female / 2, nb_marker_male / 2)
    print "%s discarded because not supported by the format (%s)" % (non_supported_markers,
                                                                     ', '.join(["%s:%s" % (genotype, count) for
                                                                                genotype, count in
                                                                                non_supported_genotypes.items()]))
    print "%s discarded because of bad pattern in offspring" % bad_pattern


def output_marker_in_joinmap(name, possible_genotypes, ordered_samples, geno_parents, all_offsprings_genotypes,
                             max_missing_offspring):
    possible_offspring_genotypes = possible_genotypes.get(geno_parents)
    if not possible_offspring_genotypes:
        return None
    all_sample_out = []
    bad_offsprings = 0

    for i, sample in enumerate(ordered_samples):
        geno = all_offsprings_genotypes[i]
        result = possible_offspring_genotypes.get(geno)
        if result:
            all_sample_out.append(result)
        else:
            all_sample_out.append("--")
    if "".join(all_sample_out).count("-") / 2 > max_missing_offspring:
        logging.info(
            "nb missing offspring = %s>%s: %s" % ("".join(all_sample_out).count("-") / 2, max_missing_offspring, name))
        return ""
    else:
        return "%s\t%s\t%s\n" % (name, possible_genotypes.get(geno_parents).get('type'), ' '.join(all_sample_out))


def output_marker_in_onemap(name, possible_genotypes, ordered_samples, geno_parents, all_offsprings_genotypes,
                            max_missing_offspring):
    possible_offspring_genotypes = possible_genotypes.get(geno_parents)
    if not possible_offspring_genotypes:
        return
    all_sample_out = []

    for i, sample in enumerate(ordered_samples):
        geno = all_offsprings_genotypes[i]
        result = possible_offspring_genotypes.get(geno)
        if result:
            all_sample_out.append(result)
        else:
            all_sample_out.append("--")
    if "".join(all_sample_out).count("-") / 2 > max_missing_offspring:
        logging.info(
            "nb missing offspring = %s>%s: %s" % ("".join(all_sample_out).count("-") / 2, max_missing_offspring, name))
        return ""
    else:
        return "*%s\t%s\t%s\n" % (name, possible_genotypes.get(geno_parents).get('type'), ','.join(all_sample_out))


def output_marker_in_cartagene(name, possible_genotypes, ordered_samples, geno_parents, all_offsprings_genotypes):
    female_genotypes, male_genotypes = possible_genotypes
    # female map
    out_female = []
    #female_out.write("*%s\t%s\n"%(marker_name, key1))
    possible_offspring_genotypes = female_genotypes.get(geno_parents)
    if not possible_offspring_genotypes:
        out_female = None
    else:
        for i, sample in enumerate(ordered_samples):
            geno = all_offsprings_genotypes[i]
            result = possible_offspring_genotypes.get(geno)
            if result:
                out_female.append(result)
            else:
                out_female.append("-")
    #female map 
    out_male = []
    #female_out.write("*%s\t%s\n"%(marker_name, key1))
    possible_offspring_genotypes = male_genotypes.get(geno_parents)
    if not possible_offspring_genotypes:
        out_male = None
    else:
        for i, sample in enumerate(ordered_samples):
            geno = all_offsprings_genotypes[i]
            result = possible_offspring_genotypes.get(geno)
            if result:
                out_male.append(result)
            else:
                out_male.append("-")
    if out_female is None:
        return_female = None
    elif "".join(out_female).count("-") > len(out_female) - 4:
        return_female = ""
    else:
        return_female = "*%s\t%s\n*%s_rev\t%s\n" % (name, "".join(out_female), name, reverse_consensus(out_female))

    if out_male is None:
        return_male = None
    elif "".join(out_male).count("-") > len(out_male) - 4:
        return_male = ""
    else:
        return_male = "*%s\t%s\n*%s_rev\t%s\n" % (name, "".join(out_male), name, reverse_consensus(out_male))

    return return_female, return_male


def compare_offspring_haplotype_to_parents(phased_genotype_mother, phased_genotype_father,
                                           phased_vcf_record, samples,
                                           genotype_quality_threshold=20,
                                           max_missing_offspring=5,
                                           perfect_match=False,
                                           allow_mendelian_error=True):
    """This Method will compare the offspring haplotype to the parental haplotype.
    It returns the single digit haplotype code genotypes that correspond to one or several SNP.
    If allow_mendelian_error is set then it will report a separate code for things that looks like mendelian errors."""
    # get the parents haplotypes
    set_parents_haplotypes = Counter()
    hm1, hm2 = merge_phased_genotype(phased_genotype_mother)
    set_parents_haplotypes[hm1] += 1
    set_parents_haplotypes[hm2] += 1

    hf1, hf2 = merge_phased_genotype(phased_genotype_father)
    set_parents_haplotypes[hf1] += 1
    set_parents_haplotypes[hf2] += 1
    if len(set_parents_haplotypes) < 2:
        logging.debug("No variant in Parental genotype: abort matching")
        return '%s %s' % (phased_genotype_mother, phased_genotype_father), None
    #sort the haplotype by the most frequent in the parent so most frequent is called 0 and second most frequent is 1 and so on
    sorted_set = sorted(set_parents_haplotypes.iteritems(), key=operator.itemgetter(1), reverse=True)
    list_parents_haplotypes = [key for key, value in sorted_set]
    parents_haplotypes = {}
    i = 0
    for hap in list_parents_haplotypes:
        parents_haplotypes[hap] = i
        i += 1
    mendelian_error_allele = i
    #Set the parental genotype based on the defined haplotypes
    genotype_mother = [str(parents_haplotypes.get(hm1)), str(parents_haplotypes.get(hm2))]
    genotype_mother.sort()
    genotype_father = [str(parents_haplotypes.get(hf1)), str(parents_haplotypes.get(hf2))]
    genotype_father.sort()
    all_offspring_genotype_high_qual = get_genotype_for_list_of_SNPs(list_vcf_record=phased_vcf_record,
                                                                     samples=samples,
                                                                     genotype_quality_threshold=genotype_quality_threshold)
    missing_offspring = 0
    all_samples_genotypes = []
    #print parents_haplotypes
    logging.debug("Set parental haplotypes to " + str(parents_haplotypes) + " where mother is " + str(
        genotype_mother) + " and father is " + str(genotype_father))

    for i in range(len(samples)):
        non_informative = False
        haplo_off1, haplo_off2 = merge_phased_genotype(all_offspring_genotype_high_qual[i])
        allele1 = parents_haplotypes.get(haplo_off1)
        allele2 = parents_haplotypes.get(haplo_off2)

        #Check for the presence of confident haplotype that shows Mendelian error
        m_haplo_off1, m_haplo_off2 = match_haplotypes_to_parents(list_parents_haplotypes, haplo_off1, haplo_off2)
        if m_haplo_off1 is None or m_haplo_off2 is None:
            if m_haplo_off1 is None:
                logging.debug('Mendelian error in haplotype %s' % haplo_off1)
                if allow_mendelian_error:
                    allele1 = mendelian_error_allele
            if m_haplo_off2 is None:
                logging.debug('Mendelian error in haplotype %s' % haplo_off2)
                if allow_mendelian_error:
                    allele2 = mendelian_error_allele

        if haplo_off1.count('.') == len(haplo_off1) or haplo_off2.count('.') == len(haplo_off2):
            logging.debug("For " + samples[
                i] + " confident haplotypes: " + haplo_off1 + " " + haplo_off2 + " have no information")

        else:
            m_haplo_off1 = m_haplo_off2 = '-' * len(haplo_off1)
            logging.debug("For " + samples[
                i] + " confident haplotypes: " + haplo_off1 + " " + haplo_off2 + " with perfect match to " + str(
                allele1) + "-" + str(allele2))
            if (allele1 is None or allele2 is None) and not perfect_match:
                #Check for imperfect matches
                m_haplo_off1, m_haplo_off2 = match_haplotypes_to_parents(list_parents_haplotypes, haplo_off1,
                                                                         haplo_off2)
                if allele1 is None: allele1 = parents_haplotypes.get(m_haplo_off1)
                if allele2 is None: allele2 = parents_haplotypes.get(m_haplo_off2)
                logging.debug("For " + samples[
                    i] + " confident haplotypes: " + haplo_off1 + " " + haplo_off2 + " with imperfect match to " + str(
                    allele1) + "-" + str(allele2))
            #This part assume mendelian segregation of the markers and try to find a match using information of one of the allele
            if (allele1 is None or allele2 is None) and not perfect_match:
                if allele1 != None:
                    if str(allele1) in genotype_mother and str(allele1) in genotype_father:
                        remaining_haplo = []
                        for h in parents_haplotypes.keys():
                            if h != m_haplo_off1:
                                remaining_haplo.append(h)
                        if len(remaining_haplo) > 0:
                            dummy, m_haplo_off2 = match_haplotypes_to_parents(remaining_haplo, '', haplo_off2)
                        allele2 = parents_haplotypes.get(m_haplo_off2)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele1) + " is both maternal and paternal, imperfect match " + haplo_off2 + " to " + str(
                            remaining_haplo) + " resolve to " + str(allele2))
                    elif str(allele1) in genotype_mother:
                        dummy, m_haplo_off2 = match_haplotypes_to_parents([hf1, hf2], '', haplo_off2)
                        allele2 = parents_haplotypes.get(m_haplo_off2)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele1) + " is maternal, imperfect match " + haplo_off2 + " to [" + hf1 + "," + hf2 + "] resolve to " + str(
                            allele2))
                    elif str(allele1) in genotype_father:
                        dummy, m_haplo_off2 = match_haplotypes_to_parents([hm1, hm2], '', haplo_off2)
                        allele2 = parents_haplotypes.get(m_haplo_off2)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele1) + " is paternal, imperfect match " + haplo_off2 + " to [" + hm1 + "," + hm2 + "] resolve to " + str(
                            allele2))
                elif allele2 != None:
                    if str(allele2) in genotype_mother and str(allele2) in genotype_father:
                        remaining_haplo = []
                        for h in parents_haplotypes.keys():
                            if h != m_haplo_off2:
                                remaining_haplo.append(h)
                        if len(remaining_haplo) > 0:
                            m_haplo_off1, dummy = match_haplotypes_to_parents(remaining_haplo, haplo_off1, '')
                        allele1 = parents_haplotypes.get(m_haplo_off1)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele2) + " is both maternal and paternal, imperfect match " + haplo_off1 + " to " + str(
                            remaining_haplo) + " resolve to " + str(allele1))
                    elif str(allele2) in genotype_mother:
                        m_haplo_off1, dummy = match_haplotypes_to_parents([hf1, hf2], haplo_off1, '')
                        allele1 = parents_haplotypes.get(m_haplo_off1)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele2) + " is maternal, imperfect match " + haplo_off1 + " to [" + hf1 + "," + hf2 + "] resolve to " + str(
                            allele1))
                    elif str(allele2) in genotype_father:
                        m_haplo_off1, dummy = match_haplotypes_to_parents([hm1, hm2], haplo_off1, '')
                        allele1 = parents_haplotypes.get(m_haplo_off1)
                        logging.debug("For " + samples[i] + " Knowing " + str(
                            allele2) + " is paternal, imperfect match " + haplo_off1 + " to [" + hm1 + "," + hm2 + "] resolve to " + str(
                            allele1))
            #This also assume mendelian segregation and check in independent parents the presence of the same allele 
            if (allele1 is None and allele2 is None) and not perfect_match and haplo_off1 == haplo_off2:
                m_haplo_off1, dummy = match_haplotypes_to_parents([hm1, hm2], haplo_off1, '')
                m_haplo_off2, dummy = match_haplotypes_to_parents([hf1, hf2], haplo_off1, '')
                allele1 = parents_haplotypes.get(m_haplo_off1)
                allele2 = parents_haplotypes.get(m_haplo_off2)
                logging.debug("For " + samples[
                    i] + " Because " + haplo_off1 + "=" + haplo_off2 + " mathing them to independent parent resolve to " + str(
                    allele1) + " and " + str(allele2))

        if allele1 is None or allele2 is None:
            missing_offspring += 1
            if allele1 is None: allele1 = '.'
            if allele2 is None: allele2 = '.'
            logging.debug("For " + samples[i] + " matching failed")

        genotype_offspring = [str(allele1), str(allele2)]
        genotype_offspring.sort()
        all_samples_genotypes.append('/'.join(genotype_offspring))
    geno_parents = '%s %s' % ('/'.join(genotype_mother), '/'.join(genotype_father))
    logging.debug("nb missing offspring = %s" % (missing_offspring))
    if missing_offspring > max_missing_offspring:
        all_samples_genotypes = None
        name = get_list_of_SNPs_name(phased_vcf_record)
        logging.info("nb missing offspring = %s>%s: %s" % (missing_offspring, max_missing_offspring, name))
    #print all_samples_genotypes
    return geno_parents, all_samples_genotypes


def match_haplotype_to_parents(list_parents_haplotypes, haplotype_offspring):
    """This method take one haplotype from an offspring and match them against the parental haplotypes.
    @return the matching parental haplotypes, the offspring if the parental haplotype wasn't found or None if Mendelian error"""
    all_results = {}
    # print 'parent haplotypes = %s'%(' '.join(list_parents_haplotypes))
    for haplotype in list_parents_haplotypes:
        res = hamming1(haplotype, haplotype_offspring)
        if all_results.has_key(res):
            all_results[res].append(haplotype)
        else:
            all_results[res] = [haplotype]

    scores = all_results.keys()
    scores.sort(reverse=False)
    #print '%s --> %s best score =%s (%s)'%(haplotype_offspring,len(all_results.get(scores[0])), scores[0], ', '.join([str(val) for val in scores]))
    if len(all_results.get(scores[0])) > 1:
        #Can't discriminate between two genotypes
        return haplotype_offspring
    elif scores[0] >= mendelian_error:
        #Mendelian error
        return None
    else:
        return all_results.get(scores[0])[0]


def match_haplotypes_to_parents(list_parents_haplotypes, ho1, ho2):
    """This method take two haplotypes from an offspring and match them against the parental haplotypes.
    @return the matching parental haplotypes or the offspring if the parental haplotype wasn't found."""
    ret1 = match_haplotype_to_parents(list_parents_haplotypes, ho1)
    ret2 = match_haplotype_to_parents(list_parents_haplotypes, ho2)
    return ret1, ret2


def hamming1(str1, str2):
    return sum(itertools.imap(not_equal, str1, str2))
    # return sum(itertools.imap(str.__ne__, str1, str2))


def not_equal(str1, str2):
    """This function compare two single position haplotypes and return 0 in case of equality 1 in case of missing genotype and 50 in case of conflicting haplotype"""
    if str1 == '.' or str2 == '.':
        return 1
    elif str1 == str2:
        return 0
    else:
        return mendelian_error


def get_genotype_for_list_of_SNPs(list_vcf_record, samples, genotype_quality_threshold=0, depth_threshold=0):
    """This script will create a list of genotype one for each offsping removing any genotype that does not pass the genotype quality of depth thresholds."""
    all_samples_genotype = []
    for sample in samples:
        all_samples_genotype.append([])

    for vcf_record in list_vcf_record:
        genotypes = vcf_record.get_all_genotype(sample_list=samples)
        genotype_qualities = vcf_record.get_all_genotype_quality(sample_list=samples)
        genotype_depth = vcf_record.get_all_sample_depth(sample_list=samples)
        for i in range(len(samples)):
            if genotype_qualities[i] is not None and float(genotype_qualities[i]) > genotype_quality_threshold and \
                            int(genotype_depth[i]) > depth_threshold:
                all_samples_genotype[i].append(re.split(r'[|/]', genotypes[i]))
            else:
                all_samples_genotype[i].append(['.', '.'])
    return all_samples_genotype


def get_list_of_SNPs_name(list_vcf_record):
    ref = None
    positions = []
    for vcf_record in list_vcf_record:
        if ref and ref != vcf_record.get_reference():
            logging.error("%s and %s are different: can only get a name to SNPs of the same reference" % (
            ref, vcf_record.get_reference()))
            return None
        else:
            ref = vcf_record.get_reference()
        positions.append(str(vcf_record.get_position()))
    return '%s:%s' % (ref, '_'.join(positions))


def merge_phased_genotype(phased_genotype):
    first = []
    second = []
    for f, s in phased_genotype:
        first.append(f)
        second.append(s)
    return ''.join(first), ''.join(second)


def add_to_array_in_dictionnary(dict, key, value):
    array = dict.get(key)
    if not array:
        array = []
        dict[key] = array
    array.append(value)


def snps_to_pattern(vcf_file, sex_info, mother, father, geno_qual_threshold):
    sample_to_sex, sex_to_sample, ordered_sample = load_sex_info(sex_info)
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader = vcfIO.VcfReader(file_handle)
    sample_names = reader.get_sample_names()
    discard_quality_parent = 0
    for vcf_records in reader:
        # First check that the parent are callable
        geno_mother = vcf_records.get_genotype(mother)
        gq_mother = vcf_records.get_genotype_quality(mother)
        geno_father = vcf_records.get_genotype(father)
        gq_father = vcf_records.get_genotype_quality(father)
        if gq_mother < geno_qual_threshold or gq_father < geno_qual_threshold:
            discard_quality_parent += 1
            continue
        #print "%s %s"%(geno_mother,geno_father)
        all_sample_out = []
        sum = 0
        if "%s %s" % (geno_mother, geno_father) == "0/1 0/1":
            continue
        for male_sample in sex_to_sample.get("M"):
            g = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=10,
                                                          minimum_depth=6, sample_list=[male_sample])
            if len(g) > 0:
                geno = g.keys()[0]
                if geno == geno_father:
                    all_sample_out.append("1")
                    sum += 1
                else:
                    all_sample_out.append("0")
            else:
                all_sample_out.append("-")
        for female_sample in sex_to_sample.get("F"):
            g = vcf_records.get_valid_genotype_per_sample(genotype_quality_threshold=10,
                                                          minimum_depth=6, sample_list=[female_sample])
            if len(g) > 0:
                geno = g.keys()[0]
                if geno == geno_mother:
                    all_sample_out.append("1")
                    sum += 1
                else:
                    all_sample_out.append("0")
            else:
                all_sample_out.append("-")

        print vcf_records.get_reference(), "".join(all_sample_out) + " %s" % sum


def filter_alleles(sample_to_allele):
    all_alleles_count = sample_to_allele.pop('all')
    total = 0.0
    valid_alleles = []
    for alleles in all_alleles_count:
        total += all_alleles_count.get(alleles)
    for allele in all_alleles_count:
        if all_alleles_count.get(allele) > 1 and all_alleles_count.get(allele) / total > 0.05:
            valid_alleles.append(allele)
    for sample in sample_to_allele.keys():
        test_alleles = sample_to_allele.get(sample)
        total = 0.0
        for alleles in test_alleles:
            total += test_alleles.get(alleles)
        for allele in test_alleles.keys():
            if not allele in valid_alleles or test_alleles.get(allele) == 1 or test_alleles.get(allele) / total <= 0.05:
                test_alleles.pop(allele)


def detect_missing_allele(vcf_file, normalization_factors):
    file_handle = utils_logging.open_input_file(vcf_file, pipe=False)
    reader = vcfIO.VcfReader(file_handle)
    reader.get_header_lines()
    mother = "MP"
    father = "FP"
    for vcf_recors in reader:
        d_mo = vcf_recors.get_sample_depth(mother)
        d_fa = vcf_recors.get_sample_depth(father)
        norm_mother = float(d_mo) / normalization_factors.get(mother)
        norm_father = float(d_fa) / normalization_factors.get(father)
        print norm_mother, norm_father, vcf_recors.get_reference(), vcf_recors.get_position()


def main():
    # initialize the logging
    utils_logging.init_logging()
    #Setup options
    optparser = _prepare_optparser()
    (options, args) = optparser.parse_args()
    #verify options
    arg_pass = _verifyOption(options)
    if not arg_pass:
        logging.warning(optparser.get_usage())
        logging.critical("Non valid arguments: exit")
        sys.exit(1)
    if options.debug:
        utils_logging.init_logging(logging.DEBUG)
    if options.type_output == OUTPUT_TYPE_JOINMAP:
        phased_snps_to_joinmap(options.input_vcf_file, options.sex_info_file, options.mother_name,
                               options.father_name, geno_qual_threshold=options.geno_qual_threshold,
                               output_markers=options.output_marker, output_vcf=None,
                               max_missing_offspring=options.max_missing_offspring)
    elif options.type_output == OUTPUT_TYPE_ONEMAP:
        phased_snps_to_onemap(options.input_vcf_file, options.sex_info_file, options.mother_name,
                              options.father_name, geno_qual_threshold=options.geno_qual_threshold,
                              output_markers=options.output_marker, output_vcf=None,
                              max_missing_offspring=options.max_missing_offspring)
    elif options.type_output == OUTPUT_TYPE_CARTHAGENE:
        phased_snps_to_carthagene(options.input_vcf_file, options.sex_info_file, options.mother_name,
                                  options.father_name, geno_qual_threshold=options.geno_qual_threshold,
                                  output_markers_female=options.output_female_marker,
                                  output_markers_male=options.output_male_marker,
                                  output_vcf=options.output_vcf, max_missing_offspring=options.max_missing_offspring)


def _prepare_optparser():
    """Prepare optparser object. New options will be added in this
    function first.
    """
    usage = """usage: %prog -i <input_vcf_file> -s <sex_info_file> -m <mother_name> -f <father_name> -o <output_markers> [-g 20]"""
    description = """This script ."""

    optparser = OptionParser(description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h", "--help", action="help", help="show this help message and exit.")
    optparser.add_option("-i", "--input_vcf_file", dest="input_vcf_file", type="string",
                         help="Path to input vcf file where the SNPs are located. Default: %default")
    optparser.add_option("-s", "--sex_info_file", dest="sex_info_file", type="string",
                         help="""Path to sex information file where the name and sex of the samples is specified in tab separated format.
The format is: name<tab>M/F.
Default: %default""")
    optparser.add_option("-m", "--mother_name", dest="mother_name", type="string",
                         help="The name of the mother. Default: %default")
    optparser.add_option("-f", "--father_name", dest="father_name", type="string",
                         help="The name of the father. Default: %default")
    optparser.add_option("-g", "--geno_qual_threshold", dest="geno_qual_threshold", type="int", default=20,
                         help="the genotype quality threshold above which genotypes will be used. Default: %default")
    optparser.add_option("-o", "--output_marker", dest="output_marker", type="string",
                         help="Path to the file that will contain the markers. Default: %default")
    optparser.add_option("-p", "--output_male_marker", dest="output_male_marker", type="string",
                         help="Path to the file that will contain the males markers. Default: %default")
    optparser.add_option("-q", "--output_female_marker", dest="output_female_marker", type="string",
                         help="Path to the file that will contain the female markers. Default: %default")
    optparser.add_option("-v", "--output_vcf", dest="output_vcf", type="string",
                         help="Path to the file that will contain the SNPs retained in VCF format. Default: %default")
    optparser.add_option("-x", "--max_missing_offspring", dest="max_missing_offspring", type="int",
                         help="the maximum of missing offspring allow before a marker get removed. Offspring can be missing because the genotype is imprecise of because of a mendelian error. Default: %default")
    optparser.add_option("-t", "--type_output", dest="type_output", type="string", default=OUTPUT_TYPE_JOINMAP,
                         help="Type of output format required. (Should be one of " + ", ".join(
                             OUTPUT_TYPE) + ") Default: %default")
    optparser.add_option("--debug", dest="debug", action="store_true", default=False,
                         help="Set the verbosity ot debug mode. Default: %default")
    return optparser


def _verifyOption(options):
    """Check if the mandatory option are present in the options objects.
    @return False if any argument is wrong."""
    arg_pass = True

    if not options.input_vcf_file or not os.path.exists(options.input_vcf_file):
        logging.error("You must specify a valid input file.")
        arg_pass = False
    if not options.sex_info_file or not os.path.exists(options.sex_info_file):
        logging.error("You must specify a valid sex information file.")
        arg_pass = False
    if not options.mother_name:
        logging.error("You must specify the name of the mother.")
        arg_pass = False
    if not options.father_name:
        logging.error("You must specify the name of the father.")
        arg_pass = False
    if not options.type_output in OUTPUT_TYPE:
        logging.error("You must specify a valid output type (%s)." % ", ".join(OUTPUT_TYPE))
    elif (
            options.type_output == OUTPUT_TYPE_JOINMAP or options.type_output == OUTPUT_TYPE_ONEMAP) and not options.output_marker:
        logging.error("You must specify an output file with -o for joinmap format output.")
        arg_pass = False
    elif options.type_output == OUTPUT_TYPE_CARTHAGENE:
        if not options.output_male_marker:
            logging.error("You must specify a male output file with -p for carthagene format output.")
            arg_pass = False
        if not options.output_male_marker:
            logging.error("You must specify a female output file with -q for carthagene format output.")
            arg_pass = False
    return arg_pass


if __name__ == "__main__":
    main()

if __name__ == "1__main__":
    vcf_file = sys.argv[1]
    normalization_factors = {"MP": 1, "FP": 1.1194}
    detect_missing_allele(vcf_file, normalization_factors=normalization_factors)

if __name__ == "1__main__":
    all_h = [0, 1, 2, 3]
    all_g = []
    parents = []
    for i in range(len(all_h)):
        for j in range(i, len(all_h)):
            all_g.append("%s/%s" % (all_h[i], all_h[j]))
    for g1 in all_g:
        for g2 in all_g:
            parents.append('%s %s' % (g1, g2))
    for parent in parents:
        p1, p2 = parent.split()
        haplo1 = set(p1.split('/'))
        haplo2 = set(p2.split('/'))
        haplo1.add('.')
        haplo2.add('.')
        offsprings = set()
        for h1 in haplo1:
            for h2 in haplo2:
                offsprings.add('"%s":0' % '/'.join(sorted([h1, h2])))
        print '"%s" : {%s},' % (parent, ', '.join(offsprings))
    sys.exit(0)
    list_parents_haplotypes = ['11', '10']
    ho1 = '1.'
    ho2 = '01'
    for offspring in [ho1, ho2]:
        match_haplotypes_to_parents
    res = match_haplotypes_to_parents(list_parents_haplotypes, ho1, ho2)

    print res