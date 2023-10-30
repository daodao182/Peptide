"""
Simple example script demonstrating how to use the PeptideBuilder library.

The script generates a peptide consisting of six arginines in alpha-helix
conformation, and it stores the peptide under the name "example.pdb".
"""

from PeptideBuilder import Geometry
import PeptideBuilder

# aa_helix_linker
# geo = Geometry.geometry("H")
# geo4 = Geometry.geometry("H")
# structure = PeptideBuilder.initialize_res(geo)
# PeptideBuilder.add_residue(structure, geo4)
# for i in range(3):
#     PeptideBuilder.add_residue(structure, geo)
#     PeptideBuilder.add_residue(structure, geo4)
# PeptideBuilder.add_terminal_OXT(structure)

# geo2 = Geometry.geometry("I")
# PeptideBuilder.add_residue_AA_AA(structure, geo2)
# # linker
# geo3 = Geometry.geometry("D")
# PeptideBuilder.add_residue_aa_linker2_3(structure, geo3)
#
# # secondchain
# PeptideBuilder.add_residue_linker2_3_aa(structure,geo)
# PeptideBuilder.add_residue(structure, geo4)
# for i in range(3):
#     PeptideBuilder.add_residue(structure, geo)
#     PeptideBuilder.add_residue(structure, geo4)

# aa_ala_linker
geo = Geometry.geometry("A")
structure = PeptideBuilder.initialize_res_natural(geo)
geo2 = Geometry.geometry("M")
geo3 = Geometry.geometry("O")

for i in range(15):
    PeptideBuilder.add_residue_ala_aa(structure, geo2)
    PeptideBuilder.add_residue_aa_ala(structure, geo3)
PeptideBuilder.add_terminal_OXT(structure)
# geo4 = Geometry.geometry("L")
# PeptideBuilder.add_residue_alalinker(structure, geo4)
# geo5 = Geometry.geometry("J")
# PeptideBuilder.add_residue_linker_ala(structure,geo5)
# for i in range(8):
#     PeptideBuilder.add_residue_ala_aa(structure, geo2)
#     PeptideBuilder.add_residue_aa_ala(structure, geo3)


import Bio.PDB

out = Bio.PDB.PDBIO()
out.set_structure(structure)
out.save("16ala_15aa.pdb")


