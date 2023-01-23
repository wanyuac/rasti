#!/usr/bin/env python

"""
This module defines class Sequence.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 15 Jan 2023.
"""

class Sequence:
    """ Manage a single sequence """
    def __init__(self, seqid, descr, seq):
        self.__seqid = seqid
        self.__descr = descr
        self.__seq = seq
        return
    
    @property
    def seqid(self):
        return self.__seqid
    
    @property
    def description(self):
        return self.__descr
    
    @property
    def seq(self):
        return(self.__seq)

    def set_seqid(self, new_id):
        self.__seqid = new_id
        return
