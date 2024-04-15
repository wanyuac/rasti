"""
Implementation of utility aln2mut.

This utility is a reimplementation and update of my script aln2mut.py.

Copyright (C) 2021-2024 Yu Wan <wanyuac@gmail.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 19 June 2021; the latest update: 15 Apr 2024.
"""

# Project-specific module
from module.SanityCheck import SanityCheck
from module.Alignment import Alignment

def aln2mut(aln_file, outdir, output_prefix, ref_name, vcf_to_list, var_aln):
    SanityCheck.output_dir(outdir)
    aln = Alignment(aln_file, ref_name, outdir, output_prefix)  # Initialise an Alignment object
    vcf, var_sites = aln.aln2vcf()  # Create a VCF-like table (pandas data frame) of six columns: Sample, Pos, Ref, Alt, Type (of mutation), Aux_pos
    vcf_samples = aln.vcf2mat(vcf)  # Convert a VCF data frame into a matrix of alterations (sample x variant positions) that can be aligned to a phylogenetic tree
    if vcf_to_list:
        aln.vcf2lst(vcf, vcf_samples)
    if var_aln and var_sites != None:
        aln.extract_var_sites(var_sites)
    return
