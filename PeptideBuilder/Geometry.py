"""This module is part of the PeptideBuilder library,
written by Matthew Z. Tien, Dariya K. Sydykova,
Austin G. Meyer, and Claus O. Wilke.

The Geometry module contains the default geometries of
all 20 amino acids. The main function to be used is the
geometry() function, which returns the default geometry
for the requested amino acid.

This file is provided to you under the MIT License."""

import random
from typing import List


class Geo:
    """Geometry base class"""

    residue_name: str

    # Geometry to bring together residue
    peptide_bond: float
    CA_C_N_angle: float
    C_N_CA_angle: float

    # Backbone coordinates
    N_CA_C_angle: float
    CA_N_length: float
    CA_C_length: float
    phi: float
    psi_im1: float
    omega: float

    # Carbonyl atom
    C_O_length: float
    CA_C_O_angle: float
    N_CA_C_O_diangle: float

    def __repr__(self) -> str:
        repr = ""
        for var in self.__dict__:
            repr += "%s = %s\n" % (var, self.__dict__[var])
        return repr

class AlaGeo(Geo):
    def __init__(self):

        self.CA_C_length = 1.504
        self.N_CA_C_angle = 111.85
        self.CA_N_length = 1.469

        self.C_O_length = 1.248
        self.CA_C_O_angle = 117.85
        self.N_CA_C_O_diangle = 157.11

        self.CB_CA_length = 1.517
        self.CB_CA_C_angle = 111.54
        self.CB_CA_C_N_diangle = 124.97

        self.residue_name = "A"

class AA_odd_Geo(Geo):
    def __init__(self):
        self.N_C_length = 1.293
        self.N_C_CA_angle = 119.98
        self.N_C_CA_N_diangle = -24.66

        self.CD1_N_length = 1.417
        self.CD1_N_C_angle = 123.14
        self.CD1_N_C_CA_diangle = 174.31

        self.CG_CD1_length = 1.498
        self.CG_CD1_N_angle = 108.91
        self.CG_CD1_N_C_diangle = -126.10

        self.NB_CG_length = 1.500
        self.NB_CG_CD1_angle = 112.55
        self.NB_CG_CD1_N_diangle = 65.68

        self.CA_NB_length = 1.436
        self.CA_NB_CG_angle = 114.39
        self.CA_NB_CG_CD1_diangle = 74.05

        self.C_CA_length = 1.443
        self.C_CA_NB_angle = 120.94
        self.C_CA_NB_CG_diangle = -128.72

        self.O_C_length = 1.234
        self.O_C_CA_angle = 121.09
        self.O_C_CA_NB_diangle = -176.54

        self.CE1_CD1_length = 1.493
        self.CE1_CD1_CG_angle = 109.43
        self.CE1_CD1_CG_NB_diangle = -173.67

        self.SG_NB_length = 1.634
        self.SG_NB_CG_angle = 116.07
        self.SG_NB_CG_CD1_diangle = -147.65

        self.OD1_SG_length = 1.399
        self.OD1_SG_NB_angle = 108.78
        self.OD1_SG_NB_CG_diangle = 57.86

        self.OD2_SG_length = 1.424
        self.OD2_SG_NB_angle = 105.26
        self.OD2_SG_NB_CA_diangle = -33.87

        self.CD2_SG_length = 1.733
        self.CD2_SG_NB_angle = 104.88
        self.CD2_SG_NB_CG_angle = -60.51

        self.CE2_CD2_length = 1.345
        self.CE2_CD2_SG_angle = 120.60
        self.CE2_CD2_SG_NB_diangle = -86.00

        self.CZ1_CE2_length = 1.412
        self.CZ1_CE2_CD2_angle = 119.20
        self.CZ1_CE2_CD2_SG_diangle = 176.88

        self.CE3_CD2_length = 1.418
        self.CE3_CD2_SG_angle = 120.25
        self.CE3_CD2_SG_NB_diangle = 94.20

        self.CZ2_CE3_length = 1.357
        self.CZ2_CE3_CD2_angle = 120.74
        self.CZ2_CE3_CD2_SG_diangle = -178.21

        self.CH_CZ1_length = 1.342
        self.CH_CZ1_CE2_angle = 119.53
        self.CH_CZ1_CE2_CD2_diangle = 1.24

        self.Cl17_CH_length = 1.763
        self.Cl17_CH_CZ1_angle = 119.55
        self.Cl17_CH_CZ1_CE2_diangle = 177.27

        self.residue_name = "M"

class AA_even_Geo(Geo):
    def __init__(self):
        self.N_C_length = 1.346
        self.N_C_CA_angle = 116.24
        self.N_C_CA_N_diangle = -28.08

        self.CD1_N_length = 1.487
        self.CD1_N_C_angle = 124.56
        self.CD1_N_C_CA_diangle = 173.63

        self.CG_CD1_length = 1.515
        self.CG_CD1_N_angle = 109.55
        self.CG_CD1_N_C_diangle = -118.70

        self.NB_CG_length = 1.533
        self.NB_CG_CD1_angle = 109.48
        self.NB_CG_CD1_N_diangle = 61.74

        self.CA_NB_length = 1.479
        self.CA_NB_CG_angle = 113.52
        self.CA_NB_CG_CD1_diangle = 75.63

        self.C_CA_length = 1.520
        self.C_CA_NB_angle = 116.26
        self.C_CA_NB_CG_diangle = -140.84

        self.O_C_length = 1.251
        self.O_C_CA_angle = 115.40
        self.O_C_CA_NB_diangle = -172.77

        self.CE1_CD1_length = 1.513
        self.CE1_CD1_CG_angle = 109.87
        self.CE1_CD1_CG_NB_diangle = -174.15

        self.residue_name = "N"

class Ala_odd_Geo(Geo):
    def __init__(self):

        self.N_C_length = 1.320
        self.N_C_CA_angle = 123.10
        self.N_C_CA_NB_diangle = 3.54

        self.CA_N_length = 1.421
        self.CA_N_C_angle = 121.09
        self.CA_N_C_CA_diangle = -177.47

        self.C_CA_length = 1.483
        self.C_CA_N_angle = 115.00
        self.C_CA_N_C_diangle = -72.52

        self.CB_CA_length = 1.493
        self.CB_CA_C_angle = 110.27
        self.CB_CA_C_N_diangle = 126.84


        self.C_O_length = 1.248
        self.CA_C_O_angle = 117.85
        self.N_CA_C_O_diangle = 157.19

        self.residue_name = "O"

class Ala_even_Geo(Geo):
    def __init__(self):
        self.N_C_length = 1.379
        self.N_C_CA_angle = 123.10
        self.N_C_CA_NB_diangle = 3.54

        self.CA_N_length = 1.421
        self.CA_N_C_angle = 121.09
        self.CA_N_C_CA_diangle = -177.47

        self.C_CA_length = 1.483
        self.C_CA_N_angle = 115.00
        self.C_CA_N_C_diangle = -72.52

        self.CB_CA_length = 1.493
        self.CB_CA_C_angle = 110.27
        self.CB_CA_C_N_diangle = 126.84

        self.C_O_length = 1.226
        self.CA_C_O_angle = 121.98
        self.N_CA_C_O_diangle = 152.28

        self.residue_name = "P"

#aa_ala_linker
class LinkerGeo(Geo):
    def __init__(self):
        self.N1_C_length = 1.415
        self.N1_C_CA_angle = 122.28
        self.N1_C_CA_N_diangle = -80.88

        self.C5_N1_length = 1.426
        self.C5_N1_C_angle = 129.33
        self.C5_N1_C_CA_diangle = -164.19

        self.C6_C5_length = 1.522
        self.C6_C5_N1_angle = 114.64
        self.C6_C5_N1_C_diangle = 179.90

        self.C7_C6_length = 1.539
        self.C7_C6_C5_angle = 114.64
        self.C7_C6_C5_N1_diangle = 179.05

        self.C8_C7_length = 1.480
        self.C8_C7_C6_angle = 110.54
        self.C8_C7_C6_C5_diangle = -176.02

        self.O3_C8_length = 1.224
        self.O3_C8_C7_angle = 121.04
        self.O3_C8_C7_C6_diangle = 87.98

        self.N2_C8_length = 1.336
        self.N2_C8_C7_angle = 115.83
        self.N2_C8_C7_C6_diangle = -91.43

        self.C9_N2_length = 1.447
        self.C9_N2_C8_angle = 122.53
        self.C9_N2_C8_C7_diangle = -164.62

        self.C10_C9_length = 1.497
        self.C10_C9_N2_angle = 109.17
        self.C10_C9_N2_C8_diangle = -114.51

        self.O4_C10_length = 1.227
        self.O4_C10_C9_angle = 121.22
        self.O4_C10_C9_N2_diangle = -2.91

        self.N3_C10_length = 1.341
        self.N3_C10_C9_angle = 115.99
        self.N3_C10_C9_N2_diangle = 177.40

        self.C11_N3_length = 1.435
        self.C11_N3_C10_angle = 123.13
        self.C11_N3_C10_C9_diangle = -116.00

        self.C12_C11_length = 1.527
        self.C12_C11_N3_angle = 110.53
        self.C12_C11_N3_C10_diangle = 119.92

        self.C13_C12_length = 1.535
        self.C13_C12_C11_angle = 111.39
        self.C13_C12_C11_N3_diangle = -177.32

        self.C4_C13_length = 1.502
        self.C4_C13_C12_angle = 113.91
        self.C4_C13_C12_C11_diangle = 179.83

        self.O2_C4_length = 1.212
        self.O2_C4_C13_angle = 125.03
        self.O2_C4_C13_C12_diangle = 179.27

        self.residue_name = "L"

class linker_Ala_Geo(Geo):
    def __init__(self):
        self.N_C4_length = 1.336
        self.N_C4_C13_angle = 116.74
        self.N_C4_C13_C12_diangle = 5.40

        self.CA_N_length = 1.421
        self.CA_N_C4_angle = 123.32
        # self.CA_N_C4_angle = 117.10
        self.CA_N_C4_C13_diangle = 154.74
        # self.CA_N_C4_C13_diangle = 179.92

        self.C_CA_length = 1.483
        self.C_CA_N_angle = 115.00
        self.C_CA_N_C4_diangle = -118.38
        # self.C_CA_N_C4_diangle = -70.86

        self.CB_CA_length = 1.493
        self.CB_CA_C_angle = 110.27
        self.CB_CA_C_N_diangle = 126.84

        self.C_O_length = 1.226
        self.CA_C_O_angle = 121.98
        self.N_CA_C_O_diangle = 152.28

        self.residue_name = "J"


class AAGeo(Geo):

    def __init__(self):
        self.CA_NB_length = 1.46
        self.CA_C_length = 1.52
        self.NB_CA_C_angle = 111.425

        self.CA_C_N_angle=114.536
        self.NB_CA_C_N_diangle=-146.253
        self.NB_CA_C_N_diangle1=-142.754
        self.C_N_CD1_angle=120.434

        self.C_O_length = 1.225
        self.CA_C_O_angle =119.92
        self.NB_CA_C_O_diangle1 =44.97

        self.N_CD1_length=1.427
        self.N_CD1_CG_angle=112.261
        self.N_CD1_CG_NB_diangle=66.18


        self.CD1_CG_length =1.45
        self.CD1_CG_NB_angle =114.475
        self.CD1_CG_NB_CA_diangle=-116.094
        self.CD1_CG_NB_CA_diangle1=-121.202

        self.CG_NB_length =1.54
        self.CG_NB_CA_angle =117.651
        self.CG_NB_CA_C_diangle=88.315
        self.CG_NB_CA_C_diangle1=86.201


        self.a= -121.093
        self.a1=-116.058
        self.b= 59.0
        self.c= -137.9
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.536
        self.C_N_CA_angle = 35.003


        self.NB_SG_length = 1.625
        self.CG_NB_SG_angle= 115.052
        self.CD1_CG_NB_SG_diangle=94.641
        self.CD1_CG_NB_SG_diangle1= 91.313

        self.OD2_SG_length=1.42
        self.OD2_SG_NB_angle=106.396
        #left
        self.CA_NB_SG_OD2_diangle=38.915

        self.OD1_SG_length = 1.42
        self.OD1_SG_NB_angle = 106.362
        #left
        self.CG_NB_SG_OD1_diangle=-41.980

        self.CA_C_N_CD1_diangle=-173.574
        self.CA_C_N_CD1_diangle1=-171.269

        self.CE1_CD1_length=1.530
        self.CE1_CD1_CG_angle= 108.23
        self.CE1_CD1_CG_NB_diangle = -175.64

        self.SG_CD2_length = 1.752
        self.NB_SG_CD2_angle = 108.581
        self.CG_NB_SG_CD2_diangle = 72.684

        self.CD2_CE2_length = 1.396
        self.SG_CD2_CE2_angle = 120.461
        self.NB_SG_CD2_CE2_diangle = -81.533

        self.CE2_CZ1_length = 1.392
        self.CD2_CE2_CZ1_angle = 118.296
        self.SG_CD2_CE2_CZ1_diangle = -177.738

        self.CD2_CE3_length = 1.405
        self.SG_CD2_CE3_angle = 119.029
        self.NB_SG_CD2_CE3_diangle = 99.204

        self.CE3_CZ2_length = 1.394
        self.CD2_CE3_CZ2_angle = 118.128
        self.SG_CD2_CE3_CZ2_diangle = 179.643

        self.CZ1_CH_length = 1.325
        self.CE2_CZ1_CH_angle = 120.707
        self.CD2_CE2_CZ1_CH_diangle = -1.708

        self.CH_Cl17_length = 1.750
        self.CZ1_CH_Cl17_angle = 119.309
        self.CE2_CZ1_CH_Cl17_diangle = 179.004
        #seconleft
        self.N_CL_length = 1.51
        self.N_CL_C15_angle = 106.96
        self.N_CL_C15_C14_diangle = -51.84

        self.CD1_N_CL_angle = 110.66
        self.CD1_N_CL_C15_diangle = -54.47

        self.CG_CD1_N_CL_diangle = -54.47

        #ala_aa
        self.peptide_bond_alaaa = 1.357
        self.CA_C_N_angle_alaaa = 116.89
        self.NB_CA_C_N_diangle_alaaa = -171.22

        self.C_N_CD1_angle_alaaa = 122.98
        self.CA_C_N_CD1_diangle_alaaa = 178.36

        self.c_alaaa = 79.14
        self.residue_name = "H"

class AA_AAGeo(Geo):

    def __init__(self):

        self.N_C_length = 1.51
        self.N_C_CA_angle = 120.61
        self.N_C_CA_NB_diangle = -176.69

        self.CD1_N_length = 1.516
        self.CD1_N_C_angle = 116.30
        self.CD1_N_C_CA_diangle = 111.32

        self.CG_CD1_length = 1.54
        self.CG_CD1_N_angle = 107.53
        self.CG_CD1_N_C_diangle = 171.55

        self.NB_CG_length = 1.512
        self.NB_CG_CD1_angle = 108.42
        self.NB_CG_CD1_N_diangle = 172.86

        self.CA_NB_length = 1.51
        self.CA_NB_CG_angle = 112.02
        self.CA_NB_CG_CD1_diangle = -166.18

        self.C_CA_length = 1.54
        self.C_CA_NB_angle = 112.64
        self.C_CA_NB_CG_diangle = -66.27

        self.C_O_length =1.51
        self.CA_C_O_angle = 110.40
        self.NB_CA_C_O_diangle = 63.73

        self.residue_name = "I"


def geometry(AA: str) -> Geo:
    """Generates the geometry of the requested amino acid.
    The amino acid needs to be specified by its single-letter
    code. If an invalid code is specified, the function
    returns the geometry of Glycine."""


    if AA == "H":
        return AAGeo()
    elif AA == "I":
        return AA_AAGeo()
    elif AA == "A":
        return AlaGeo()
    elif AA == "J":
        return linker_Ala_Geo()
    elif AA == "M":
        return AA_odd_Geo()
    elif AA == "N":
        return AA_even_Geo()
    elif AA == "O":
        return Ala_odd_Geo()
    elif AA == "P":
        return Ala_even_Geo()
    else:
        return LinkerGeo()
