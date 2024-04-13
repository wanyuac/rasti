#!/usr/bin/env python

"""
This module defines class Sequence.

Dependencies: Python 3

Copyright (C) 2023-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 24 Jan 2023.
"""

class Sequence:
    """ Manage a single sequence """
    def __init__(self, id, descr, seq):
        self.__id = id
        self.__descr = descr
        self.__seq = seq
        return
    
    @property
    def id(self):
        return self.__id
    
    @property
    def description(self):
        return self.__descr
    
    @property
    def seq(self):
        return(self.__seq)

    @id.setter
    def id(self, new_id):
        self.__id = new_id
        return
    
    @description.setter
    def description(self, new_descr):
        self.__descr = new_descr
        return

    @seq.setter
    def seq(self, new_seq):
        self.__seq = new_seq
        return
