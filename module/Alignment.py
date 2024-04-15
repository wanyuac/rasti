"""
Class Alignment for utility aln2mut.

Copyright (C) 2021-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 19 June 2021; the latest update: 15 Apr 2024.
"""

import os
import sys
import pandas

class Alignment:
    """ Imports an alignment and identifies mutations in this alignment """
    def __init__(self, aln_file, ref_name, outdir, output_prefix):
        self.__aln = self.__import_alignment(aln_file)  # Result: a dictionary with the structure {sequence ID : sequence}
        self.__ref_name = ref_name
        self.__outdir = outdir
        self.__output_prefix = output_prefix
        self.__ref_seq, self.__ref_seq_without_gaps, self.__aln_len = self.__validate_params(self.__aln, self.__ref_name)
        print("Info: length of the reference sequence %s: %i bp or aa" % (ref_name, len(self.__ref_seq_without_gaps)), file = sys.stderr)  # The sequence length does not count the total width of gaps.
        self.__aln = {key : value for key, value in self.__aln.items() if key != ref_name}  # Exclude the reference from the alignment, leaving an alignment of all other sequences
        return
    
    # Properties ###############
    @property
    def alignment(self):
        return self.__aln
    
    @property
    def ref_name(self):
        return self.__ref_name
    
    @property
    def outdir(self):
        return self.__outdir

    @property
    def alignment_length(self):
        return self.__aln_len

    @property
    def ref_seq(self):
        return self.__ref_seq

    @property
    def ref_seq_without_gaps(self):
        return self.__ref_seq_without_gaps
    
    @property
    def output_prefix(self):
        return self.__output_prefix

    # Public methods ###############
    def aln2vcf(self):
        """
        Mutation calling: to create a VCF-like table of five columns: Sample, Pos, Ref, Alt, Type (of mutation)
        Code for alteration types: S (substitution), I (insertion to the reference), D (deletion from the reference)
        This function implements the core algorithm of this script.
        """
        samples = list()
        coords = list()  # Positions in the reference sequence
        aux_coords = list()  # Auxiliary positions (real numbers) for sorting the output VCF file based on positions
        refs = list()  # Nucleotide bases/amino acids (AAs) in the reference sequence
        alts = list()  # Alternative bases/AAs in the sample sequence
        types = list()  # Types of alterations
        var_sites = set()  # Indices of variable sites across samples in the input alignment
        for sam, seq in self.__aln.items():  # Of note, a sample does not appear in the VCF if it is identical to the reference.
            p = 0  # A pointer for the current character in the reference sequence, excluding '-' characters. (real position - 1)
            ins = False  # A flag for whether the current position is in an insertion ('-' characters in the reference sequence)
            ins_up = 0  # Immediately upstream position of an insertion; Variables ins_up and p mark the flanking positions of the current insertion.
            ins_seq = ''  # The inserted sequence
            for i in range(0, self.__aln_len):  # Python character indexes across sequences in the alignment, including '-' characters. In the alignment file, every sequence (including '-') must have the same length.
                r = self.__ref_seq[i]  # Reference base/AA
                s = seq[i]  # Sample base/AA
                if r == s:
                    if r != '-':  # s = r = '-' when there is a larger insertion in another sequence overlapping the current sample and reference sequences.
                        p += 1
                        if ins:  # The pointer reaches the end of an insertion.
                            samples.append(sam)
                            if p == 1:
                                coords.append('^1')  # The insertion happens before the first character of the reference. Do not need to write '0^1' here.
                            else:
                                coords.append('%i^%i' % (ins_up, p))  # E.g., '25^36' means an insertion between positions 25 and 26. ins_up equals zero when the insertion happens before the start codon of the reference sequence.
                            aux_coords.append(ins_up + 0.5)  # The extra 0.5 marks the insertion site between two consecutive positions.
                            refs.append('-')
                            alts.append(ins_seq)
                            types.append('I')
                            ins = False
                            ins_seq = ''
                else:  # Record an alteration. Note that in a multisequence alignment, r and s may both be '-'.
                    var_sites.add(i)
                    if r != '-':  # The alteration is a deletion, a substitution, or the end of an insertion.
                        samples.append(sam)
                        p += 1
                        if ins:  # Now the insertion region ends
                            coords.append('^1' if p == 1 else '%i^%i' % (ins_up, p))
                            aux_coords.append(ins_up + 0.5)
                            refs.append('-')
                            alts.append(ins_seq)
                            types.append('I')
                            ins = False
                            ins_seq = ''
                        else:  # Substitution or deletion
                            coords.append(str(p))
                            aux_coords.append(p)
                            refs.append(r)
                            alts.append(s)
                            types.append('D' if s == '-' else 'S')
                    else:  # The alteration is an insertion to the reference. The variable p does not increase in this situation.
                        if ins:  # The current position remains in an insertion
                            ins_seq += s
                        else:  # At the start of an insertion: start to record the current insertion.
                            ins = True
                            ins_up = p  # Record where the insertion starts
                            ins_seq = s  # Start to record the current insertion region
            if ins:  # Finish the insertion that happens at the end of the reference sequence
                samples.append(sam)
                coords.append(str(ins_up) + '^')
                aux_coords.append(ins_up + 0.5)
                refs.append('-')
                alts.append(ins_seq)
                types.append('I')
                ins = False
                ins_seq = ''
        vcf = pandas.DataFrame({'Sample' : samples, 'Pos' : coords, 'Ref' : refs, 'Alt' : alts, 'Type' : types, 'Aux_pos' : aux_coords})
        if len(var_sites) > 0:
            var_sites = list(var_sites)
            var_sites.sort()
        else:
            print("Info: every sequence in the alignment is the same.", file = sys.stderr)
            var_sites = None
        vcf.to_csv(os.path.join(self.__outdir, self.__output_prefix + '_vcf.tsv'), index = False, sep = '\t')  # This table can be read into Python using pandas.read_csv('XXXX_vcf.tsv', sep = '\t').
        return vcf, var_sites
    
    def vcf2mat(self, vcf):
        """ Convert a VCF data frame into a matrix of alterations (sample x variant positions) """
        var_pos = vcf[['Aux_pos', 'Pos']].drop_duplicates()  # Select two columns from the data frame VCF and drop duplicated rows (There is no need to use pos.Aux_pos.unique() later as Aux_pos and Pos are linked)
        var_pos = var_pos.sort_values(by = ['Aux_pos'], ascending = True)  # Sort the data frame by auxiliary positions
        var_pos = var_pos.reset_index()
        n = len(var_pos.index)  # Row count of var_pos
        mat = pandas.DataFrame(columns = ['Sample'] + var_pos['Pos'].tolist())  # Initiate an empty matrix (as a data frame) with defined column names
        new_row = pandas.Series([self.__ref_name] + self.__get_ref_row(ref_gap_free = self.__ref_seq_without_gaps, coords = var_pos['Aux_pos'].tolist()), index = mat.columns)  # Create a new row with the reference sequence's name and genotypes/AAs at variable sites
        mat = pandas.concat([mat, new_row.to_frame().T], ignore_index = True)
        vcf_samples = vcf['Sample'].unique()  # Returns an array object
        for sam in self.__aln.keys():  # Note that keys do not include the reference sequence; add a new row into the matrix mat for each sample
            sample_row_content = [sam]  # A list of values for fields in the new row
            if sam in vcf_samples:
                vcf_sam = vcf.loc[vcf['Sample'] == sam]  # Select rows corresponding to the current sample
                vcf_sam_auxpos = vcf_sam['Aux_pos'].tolist()
                for _, row in var_pos.iterrows():  # Iterate through positions of variants (namely, non-ID columns of the output matrix)
                    p = row['Aux_pos']  # A numeric value.
                    if p in vcf_sam_auxpos:  # The current sample has an alteration at position p.
                        sample_row_content.append(vcf_sam.loc[vcf_sam['Aux_pos'] == p, 'Alt'].iloc[0])  # A single row will be extracted from vcf_sam. Then get the first (also the only one) value from the column Alt.
                    else:
                        sample_row_content.append('.')  # The current character is identical to that in the reference.
            else:
                print("Info: sequence " + sam + " is identical to the reference sequence.", file = sys.stderr)
                sample_row_content += ['.'] * n
            new_row = pandas.Series(sample_row_content, index = mat.columns)
            mat = pandas.concat([mat, new_row.to_frame().T], ignore_index = True)  # Append a sample's row to the matrix
        mat.to_csv(os.path.join(self.__outdir, self.__output_prefix + '_matrix.tsv'), index = False, sep = '\t')  # Save the variant matrix as a TSV file
        return vcf_samples
    
    def vcf2lst(self, vcf, vcf_samples):
        """ Convert a VCF file into a list of alterations and save the list as a TSV file """
        lst = pandas.DataFrame(columns = ['Sample', 'Alteration'])
        for sam in vcf_samples:
            vcf_sam = vcf.loc[vcf['Sample'] == sam]
            alt = list()
            for _, row in vcf_sam.iterrows():
                t = row['Type']
                if t == 'S':  # Substitution
                    alt.append(row['Ref'] + row['Pos'] + row['Alt'])  # E.g. V80F
                elif t == 'I':  # Insertion
                    alt.append(row['Pos'].replace('^', 'ins' + row['Alt']))  # E.g. 10insATRQ11 (The prefix 'ins' seems redundant here, but it can be used as a keyword by users for quickly identifying insertions)
                else:  # t == 'D', deletion
                    alt.append(row['Ref'] + row['Pos'] + 'del')  # E.g. R80del. Users may want to manually merge consecutive deletions into a single one.
            new_row = pandas.Series([sam, ','.join(alt)], index = lst.columns)
            lst = pandas.concat([lst, new_row.to_frame().T], ignore_index = True)  # Sample name'\t'A comma-delimited list of alterations
        lst.to_csv(os.path.join(self.__outdir, self.__output_prefix + '_list.tsv'), index = False, sep = '\t')
        return

    def extract_var_sites(self, var_sites):
        """ Remove invariable sites from the input alignment, assuming the pre-sorted list var_sites != None. """
        aln_var = dict()
        aln_var[self.__ref_name] = self.__get_polymorphisms(self.__ref_seq, var_sites)
        for sam, seq in self.__aln.items():
            aln_var[sam] = self.__get_polymorphisms(seq, var_sites)
        with open(os.path.join(self.__outdir, self.__output_prefix + '_variable_sites.aln'), 'w') as fasta:
            for sam, seq in aln_var.items():
                fasta.write('>' + sam + '\n')
                fasta.write(seq + '\n')
        return

    # Private methods ###############
    def __import_alignment(self, fasta):
        """ Import sequences from the input FASTA file in which sequences may contain '-' characters. """
        aln = dict()
        with open(fasta, 'r') as f:
            lines = f.read().splitlines()
        s = ''  # Sequence cache
        i = None  # ID of the sequence in the cache
        for line in lines:
            if line.startswith('>'):  # A new sequence record is encountered
                if i != None:
                    aln[i] = s
                    s = ''  # Reset the sequence cache
                i = line.split(' ')[0]  # Drop the sequence annotation from the header line
                i = i[1 : ]  # Drop '>' from the ID
            else:
                s += line
        aln[i] = s  # Save the last record
        return aln

    def __validate_params(self, aln, ref_name):
        """ Check whether the reference sequence is in the alignment file """
        msg = None
        # Check sequence content
        if len(aln) == 1:
            msg = "Error: Input alignment must consist of at least two sequences."
        elif not (ref_name in aln.keys()):
            msg = "Error: Reference sequence " + ref_name + " is not found in the alignment file."
        else:
            ref_seq = aln[ref_name]  # The reference sequence (may contain '-' when insertions are present in sample sequences).
            ref_gap_free = ref_seq.replace('-', '')  # Otherwise, the coordinates in function get_ref_row do not match those in the reference sequence.
            if ref_gap_free == '':
                msg = "Error: the reference sequence cannot be a gap (namely, consisting of only dash characters)."
        # Check sequence lengths
        seq_len = list()
        for s in aln.values():
            seq_len.append(len(s))
        n = max(seq_len)  # Maximum of alignment lengths
        if min(seq_len) != n:
            msg = "Error: sequences (including gaps) in the alignment must have the same length."
        if msg != None:
            print(msg, file = sys.stderr)
            sys.exit(1)
        return ref_seq, ref_gap_free, n
    
    def __get_ref_row(self, ref_gap_free, coords):
        """ A subordinate function of method vcf2mat for making a list of reference characters corresponding to alteration positions """
        ref_chars = list()
        for p in coords:
            if int(p) == p:  # Insertion: int(p) < p; for instance, p = 10.5 and int(p) = 10.
                try:
                    i = int(p) - 1
                    ref_chars += [ref_gap_free[i]]
                except IndexError:
                    print("Runtime error: index " + str(i) + " exceeds the range of the gap-free reference sequence.", file = sys.stderr)
                    sys.exit(1)
            else:
                ref_chars += ['-']
        return(ref_chars)
    
    def __get_polymorphisms(self, seq, var_sites):
        """ A subordinate function of method extract_var_sites """
        s = ''
        for i in var_sites:
            s += seq[i]
        return s
