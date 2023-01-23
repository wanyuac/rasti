#!/usr/bin/env python

"""
This module defines class Hit.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 23 Jan 2023.
"""

from collections import namedtuple
from lib.Sequence import Sequence
from Bio import SeqIO

HIT_ATTRS = ['qlen', 'slen', 'pident', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'evalue', 'bitscore']
Hit_attrs = namedtuple('Hit_attr', HIT_ATTRS)  # from 'qlen' to 'bitscore'
CODON_TABLE_11_START = ['ATG', 'ATT', 'ATC', 'ATA', 'CTG', 'GTG', 'TTG']  # Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
CODON_TABLE_11_STOP = ['TAA', 'TAG', 'TGA']  # Sometimes the last one or two bases were missing from BLAST's output when an alternative stop codon is present.

class Hit:
    """
    Store a single hit in BLAST's outputs. This function assumes that each input line follows fields defined by
    HIT_TABLE_COLS + ['sseq'] (see BLAST.py).
    """
    def __init__(self, sample, hit_line, append_sample_name):
        """
        Parse a line in BLAST's tabulated output and create a Hit object
        Parameter hit_line must not be None.
        Hit attributes: qlen (query sequence length), slen (subject sequence length), pident (percentage of identical matches), qcovhsp (query Coverage Per HSP),
        length (alignment length), mismatch (number of mismatches), gapopen (number of gap openings), qstart (start of alignment in query), qend (end of alignment in query),
        sstart (start of alignment in subject), send (end of alignment in subject), sstrand (subject Strand), evalue (expect value), bitscore (bit score).
        """
        fields = hit_line.split('\t') if hit_line != None else [None] * 17
        self.__sample = sample  # Name of the subject genome
        self.__query = fields[0]  # Name of the query sequence used for BLAST ('qseqid')
        self.__contig = fields[1]  # Name of the subject sequence (a contig in a draft genome, a complete genome, etc) ('sseqid')
        self.__hit = '@'.join([self.__query, self.__sample]) if append_sample_name else self.__query  # Hit ID. For instance, gene1@sample1.
        self.__attr = Hit_attrs(qlen = int(fields[2]), slen = int(fields[3]), pident = float(fields[4]), qcovhsp = float(fields[5]), length = int(fields[6]),\
                                mismatch = int(fields[7]), gapopen = int(fields[8]), qstart = int(fields[9]), qend = int(fields[10]), sstart = int(fields[11]),\
                                send = int(fields[12]), sstrand = '+' if fields[13] == 'plus' else '-', evalue = fields[14], bitscore = fields[15])
        if append_sample_name:
            seq_descr = '|'.join([self.__contig, str(self.__attr.sstart) + '-' + str(self.__attr.send), self.__attr.sstrand])
        else:
            seq_descr = '|'.join([self.__sample, self.__contig, str(self.__attr.sstart) + '-' + str(self.__attr.send), self.__attr.sstrand])
        self.__sseq = Sequence(seqid = self.__hit, descr = seq_descr, seq = fields[16].replace('-', ''))  # Aligned part of subject sequence, with gap characters '-' removed.
        return
    
    @property
    def hit(self):
        return self.__hit

    @property
    def sample(self):
        return self.__sample
    
    @property
    def query(self):
        return self.__query

    @property
    def contig(self):
        return self.__contig
    
    @property
    def sseq(self):
        return self.__sseq  # A Sequence object

    @property
    def attr(self):
        return self.__attr

    @property
    def attr_values(self):
        return [self.__query, self.__contig, str(self.__attr.qlen), str(self.__attr.slen), str(self.__attr.pident), str(self.__attr.qcovhsp),\
                str(self.__attr.length), str(self.__attr.mismatch), str(self.__attr.gapopen), str(self.__attr.qstart), str(self.__attr.qend),\
                str(self.__attr.sstart), str(self.__attr.send), self.__attr.sstrand, self.__attr.evalue, self.__attr.bitscore]

    def set_hitid(self, new_id):
        self.__hit = new_id
        self.__sseq.set_seqid(new_id = self.__hit)
        return

    def write_seq(self, fasta_handle):
        """
        Write sequences into an opened file without wrapping.
        Method write_record of FastaIO.FastaWriter does not work when multiple file handles are kept open.
        """
        print('>%s %s' % (self.__sseq.seqid, self.__sseq.description), file = fasta_handle)
        print(self.__sseq.seq, file = fasta_handle)
        return
    
    def extend_cds(self, q, f):
        """
        Extend the current CDS to recover alternative start and stop codons
        Parameters: q (name of the query sequence), f (path to the subject FASTA file)
        """
        if self.__attr.qcovhsp < 100:
            qstart = self.__attr.qstart
            qend = self.__attr.qend
            qlen = self.__attr.qlen
            start_ext = (qstart > 1 and qstart <= 4)
            end_ext = (qend < qlen and qend >= (qlen - 2))
            if start_ext or end_ext:
                contig = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))[self.__contig]
                if start_ext:
                    self.__sseq = self.__extend_start(sstart = self.__attr.sstart, send = self.__attr.send, strand = self.__attr.sstrand,\
                                                      contig = contig, by = qstart - 1)
                if end_ext:
                    self.__sseq = self.__extend_end(sstart = self.__attr.sstart, send = self.__attr.send, strand = self.__attr.sstrand,\
                                                    contig = contig, by = qlen - qend)
        return

    def __extend_start(self, sstart, send, strand, contig, by):
        """
        A subordinate function for method 'extend_cds'
        This function extends the hit towards upstream and returns a Sequence object.
        Parameters:
          sstart, original start position in the subject sequence
          send, original end position in the subject sequence
          contig: subject contig sequence
          by, number of bases to extend
        """
        return

    def __extend_end(self, sstart, send, strand, contig, by):
        """
        A subordinate function for method 'extend_cds'
        This function extends the hit towards downstream and returns a Sequence object.
        Parameters:
          sstart, original start position in the subject sequence
          send, original end position in the subject sequence
          contig: subject contig sequence
          by, number of bases to extend
        """
        return
