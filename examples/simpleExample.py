"""
Simple example script demonstrating how to use the PeptideBuilder library.

The script generates a peptide consisting of six arginines in alpha-helix
conformation, and it stores the peptide under the name "example.pdb".
"""

from PeptideBuilder import Geometry
import PeptideBuilder
#aa_helix
geo=Geometry.geometry("H")
structure = PeptideBuilder.initialize_res(geo)
for i in range(7):
    PeptideBuilder.add_residue(structure,geo)
# geo = Geometry.geometry("C")
#
# structure = PeptideBuilder.initialize_res(geo)
#
# # for i in range(10):
# PeptideBuilder.add_residue(structure, geo)
# for i in range(3):
#     PeptideBuilder.add_residue1(structure, geo)
#     PeptideBuilder.add_residue(structure,geo)
# #linker
# geo2=Geometry.geometry("J")
# PeptideBuilder.add_residue2(structure,geo2)
# geo3=Geometry.geometry("C")
# PeptideBuilder.add_residue3(structure,geo3)
# #secondchain
# PeptideBuilder.add_residue4(structure,geo)
# PeptideBuilder.add_residue(structure, geo)
# for i in range(3):
#     PeptideBuilder.add_residue1(structure, geo)
#     PeptideBuilder.add_residue(structure,geo)
#
# # add terminal oxygen (OXT) to the final glycine
# PeptideBuilder.add_terminal_OXT(structure)

#aa_gly_aa_gly
geo=Geometry.geometry("H")
structure = PeptideBuilder.initialize_res(geo)
geo2=Geometry.geometry("G")
PeptideBuilder.add_residue_aagly(structure,geo2)
for i in range(3):
    PeptideBuilder.add_residue_glyaa(structure, geo)
    PeptideBuilder.add_residue_aagly(structure, geo2)

import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("aa_gly.pdb")
