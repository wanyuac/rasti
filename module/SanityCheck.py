#!/usr/bin/env python
"""
Utility class SanityCheck for checking runtime settings

Copyright (C) 2023-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 14 Jan 2023; the latest update: 15 Apr 2024.
"""

import os
import sys

class SanityCheck:
    """ Utility functions for checking settings """
    
    @staticmethod
    def output_dir(d):
        """ Creating a directory if it does not exist. """
        if os.path.exists(d):
            print(f"Directory {d} exists.", file = sys.stdout)
        else:
            print(f"Create directory {d}", file = sys.stdout)
            os.makedirs(d)
        return

    @staticmethod
    def output_file(f, message = True):
        """ Check if a file exists """
        e = os.path.exists(f)
        if (not e) and message:
            print(f"Error: file {f} does not exists.", file = sys.stderr)
        return e

    @staticmethod
    def fasta_files(fs, suf = '.fna'):
        """ Check existance of assemblies and extract sample names from filenames """
        samples = dict()  # {sample name : path to its FASTA file}
        for f in fs:
            if os.path.exists(f):
                i = os.path.basename(f).replace(suf, '')
                samples[i] = f
            else:
                print(f"Warning: Subject assembly {f} does not exist and will not be used for sequence search.", file = sys.stderr)
        return samples

    @staticmethod
    def parameter_range(v, v_min = 70.0, v_max = 100.0, v_reset = 70.0, n = 'min_identity'):
        """ Check if minimum <= v <= maximum """
        if v < v_min:
            print(f"Warning: {n} cannot be less than {v_min}. Reset its value to {v_reset}", file = sys.stderr)
            v = v_reset
        elif v > v_max:
            print(f"Warning: {n} cannot be greater than {v_max}. Reset its value to {v_reset}", file = sys.stderr)
            v = v_reset
        return v
