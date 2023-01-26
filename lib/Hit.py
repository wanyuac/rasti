#!/usr/bin/env python

"""
This module defines class Hit.

Dependencies: Python 3

Copyright (C) 2023 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 15 Jan 2023; the latest update: 26 Jan 2023.
"""

from collections import namedtuple
from lib.Sequence import Sequence
from Bio import SeqIO
import pandas as pd

HIT_ATTRS = ['qlen', 'slen', 'pident', 'qcovhsp', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'evalue', 'bitscore']  # hslen: length of the hit in the subject sequence (without gaps)
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
        self.__id = '@'.join([self.__query, self.__sample]) if append_sample_name else self.__query  # Hit ID. For instance, gene1@sample1.
        self.__attr = pd.DataFrame(columns = HIT_ATTRS.append('hslen'))  # A single-row data frame with column names starting from 'qlen' to 'bitscore' in HIT_ATTRS
        if fields[13] == 'plus':
            sstrand = '+'
            sstart = int(fields[11])
            send = int(fields[12])
        else:
            sstrand = '-'
            sstart = int(fields[12])  # This is the only difference between coordinates used in hit tables and those in raw BLAST outputs.
            send = int(fields[11])
        hslen = send - sstart + 1  # This value may be smaller than 'length' when there are gaps in the aligned subject sequence.
        self.__attr.loc[0] = [int(fields[2]), int(fields[3]), float(fields[4]), float(fields[5]), int(fields[6]), int(fields[7]), int(fields[8]),\
                              int(fields[9]), int(fields[10]), sstart, send, sstrand, hslen, fields[14], fields[15]]
        if append_sample_name:
            seq_descr = '|'.join([self.__contig, str(self.__attr['sstart'].iloc[0]) + '-' + str(self.__attr['send'].iloc[0]), str(self.__attr['sstrand'].iloc[0]), f'{hslen}bp'])
        else:
            seq_descr = '|'.join([self.__sample, self.__contig, str(self.__attr['sstart'].iloc[0]) + '-' + str(self.__attr['send'].iloc[0]), str(self.__attr['sstrand'].iloc[0]), f'{hslen}bp'])
        self.__sseq = Sequence(id = self.__id, descr = seq_descr, seq = fields[16].replace('-', ''))  # Aligned part of subject sequence, with gap characters '-' removed.
        self.__extended = False  # A flag indicating whether the hit has been extended.
        self.__append_sample_name = append_sample_name
        return
    
    @property
    def id(self):
        return self.__id

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
    def extended(self):
        return self.__extended  # True or False

    @id.setter
    def id(self, new_id):
        self.__id = new_id
        self.__sseq.id = new_id
        return
    
    def attr_values(self, extended):
        attrs = list(self.__attr.loc[0])
        if extended:
            attrs = attrs[0 : (len(attrs) - 2)]  # E-value and bit scores are not valid for the extended hit.
        return list(map(str, attrs))  # Apply str() to every element in attrs

    def write_seq(self, fasta_handle):
        """
        Write sequences into an opened file without wrapping.
        Method write_record of FastaIO.FastaWriter does not work when multiple file handles are kept open.
        """
        print('>%s %s' % (self.__sseq.id, self.__sseq.description), file = fasta_handle)
        print(self.__sseq.seq, file = fasta_handle)
        return
    
    def extend_cds(self, f):
        """
        Extend the current CDS to recover alternative start and stop codons
        Parameter: f (path to the subject FASTA file)
        """
        if float(self.__attr['qcovhsp']) < 100:
            qstart = self.__attr['qstart'].iloc[0]
            qend = self.__attr['qend'].iloc[0]
            qlen = self.__attr['qlen'].iloc[0]
            ext_start = (qstart > 1 and qstart <= 4)
            ext_end = (qend < qlen and qend >= (qlen - 2))
            start_ext_len = 0
            end_ext_len = 0
            if ext_start or ext_end:
                contig = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))[self.__contig]  # Contig is a Bio.SeqRecord.SeqRecord object.
                sstart = self.__attr['sstart'].iloc[0]
                send = self.__attr['send'].iloc[0]
                sstrand = self.__attr['sstrand'].iloc[0]
                if ext_start:
                    sstart, send, qstart, start_ext_len = self.__extend_start(sstart = sstart, send = send, sstrand = sstrand, contig = contig, qstart = qstart)  # Extended the sequence if feasible
                if ext_end:
                    sstart, send, qend, end_ext_len = self.__extend_end(sstart = sstart, send = send, sstrand = sstrand, contig = contig, qlen = qlen, qend = qend)
                if start_ext_len > 0 or end_ext_len > 0:  # Update hit attributes
                    self.__extended = True
                    aln_len = self.__attr['length'].iloc[0] + start_ext_len + end_ext_len  # Length of the alignment overlapping the subject sequence
                    self.__attr['length'] = aln_len  # New length of the alignment against the subject sequence
                    self.__attr['qcovhsp'] = (qend - qstart + 1) / qlen * 100  # New query coverage (percent)
                    self.__attr['qstart'] = qstart
                    self.__attr['qend'] = qend
                    self.__attr['sstart'] = sstart
                    self.__attr['send'] = send
                    mismatches = self.__attr['mismatch'].iloc[0] + start_ext_len + end_ext_len
                    self.__attr['hslen'] = send - sstart + 1
                    self.__attr['mismatch'] = mismatches
                    self.__attr['pident'] = round((aln_len - mismatches - self.__attr['gapopen'].iloc[0]) / aln_len * 100, 2)
                    s = contig.seq[(sstart - 1) : send]  # s is a Seq object
                    if sstrand == '-':
                        s = s.reverse_complement()
                    self.__sseq.seq = str(s)
                    if self.__append_sample_name:
                        self.sseq.description = '|'.join([self.__contig, str(self.__attr['sstart'].iloc[0]) + '-' + str(self.__attr['send'].iloc[0]), self.__attr['sstrand'].iloc[0], f'{hslen}bp']])
                    else:
                        self.sseq.description = '|'.join([self.__sample, self.__contig, str(self.__attr['sstart'].iloc[0]) + '-' + str(self.__attr['send'].iloc[0]), self.__attr['sstrand'].iloc[0], f'{hslen}bp']])
                    ext_action = f"s{start_ext_len};t{end_ext_len}"
                else:
                    ext_action = ""
            else:
                ext_action = ""
        else:
            ext_action = ""
        return ext_action

    def __extend_start(self, sstart, send, sstrand, contig, qstart):
        """
        A private subordinate function for method 'extend_cds'
        This function extends the hit towards upstream and returns a Bio::Seq object.
        Parameters:
          sstart, original start position in the subject sequence
          send, original end position in the subject sequence
          contig: A subject contig sequence
          by, number of bases to extend
        Assumption of this function: qstart > 1.
        """
        sstart_new = sstart  # Default outputs
        send_new = send
        qstart_new = qstart
        ext_len = 0  # Number of bases (mismatches) extended
        by = qstart - 1  # Number of bases to go upstream
        if sstrand == '+':  # The first extendable condition
            i = sstart - by  # Potential new start
            if i > 0:
                codon = contig.seq[(i - 1) : (i + 2)]  # The potential start codon (a Bio.Seq.Seq object with a 3-base sequence)
                if str(codon) in CODON_TABLE_11_START:  # An alternative start codon is found; codon is a Seq object.
                    sstart_new = i
                    qstart_new = 1
                    ext_len = by
        else:
            i = send + by  # Potential new end
            if i <= len(contig.seq):  # The second extendable condition
                codon = contig.seq[(i - 3) : i]  # codon is a Bio.Seq.Seq object with a 3-base sequence
                codon = codon.reverse_complement()  # A Bio.Seq.Seq object
                if str(codon) in CODON_TABLE_11_START:
                    send = i
                    qstart_new = 1
                    ext_len = by
        return sstart_new, send_new, qstart_new, ext_len

    def __extend_end(self, sstart, send, sstrand, contig, qlen, qend):
        """
        A private subordinate function for method 'extend_cds'
        This function extends the hit towards downstream and returns a Bio::Seq object.
        Parameters:
          sstart, original start position in the subject sequence
          send, original end position in the subject sequence
          contig: subject contig sequence
          by, number of bases to extend
        Assumption of this function: qend < qlen.
        """
        sstart_new = sstart  # Default values
        send_new = send
        qend_new = qend
        ext_len = 0
        by = qlen - qend
        if sstrand == '+':
            i = send + by  # Potential new end
            if i <= len(contig.seq):
                codon = contig.seq[(i - 3) : i]  # Codon is a Seq object.
                if str(codon) in CODON_TABLE_11_STOP:
                    send_new = i
                    qend_new = qlen
                    ext_len = by
        else:
            i = sstart - by
            if i > 0:
                codon = contig.seq[(i - 1) : (i + 2)]
                codon = codon.reverse_complement()  # Codon is a Seq object.
                if str(codon) in CODON_TABLE_11_STOP:
                    sstart_new = i
                    qend_new = qlen
                    ext_len = by
        return sstart_new, send_new, qend_new, ext_len
