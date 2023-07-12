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
    """Geometry of Alanin"""

    def __init__(self):
        self.CA_N_length = 1.46
        self.CA_C_length = 1.52
        self.N_CA_C_angle = 111.068

        self.C_O_length = 1.23
        self.CA_C_O_angle = 120.5
        self.N_CA_C_O_diangle = -60.5

        self.phi = 84.22
        self.psi_im1 = 173.13
        self.omega = -170.58
        self.peptide_bond = 1.357
        self.CA_C_N_angle = 116.25
        self.C_N_CA_angle = 111.00
        
        self.CA_CB_length = 1.52
        self.C_CA_CB_angle = 109.5
        self.N_C_CA_CB_diangle = 122.6860
        
        self.residue_name = "A"


#helix_first_sidechain
class HfsGeo(Geo):
    def __init__(self):
        self.CA_NB_length = 1.46
        self.CA_C_length = 1.52
        self.NB_CA_C_angle = 111.425

        self.CA_C_N_angle = 114.536
        self.NB_CA_C_N_diangle = -146.253
        self.NB_CA_C_N_diangle1 = -142.754
        self.C_N_CD1_angle = 120.434

        self.C_O_length = 1.225
        self.CA_C_O_angle = 119.92
        self.NB_CA_C_O_diangle = 41.184
        self.NB_CA_C_O_diangle1 = 44.937

        self.N_CD1_length = 1.427
        self.N_CD1_CG_angle = 112.261
        self.N_CD1_CG_NB_diangle = 61.910
        self.N_CD1_CG_NB_diangle1 = 66.252

        self.CD1_CG_length = 1.45
        self.CD1_CG_NB_angle = 114.475
        self.CD1_CG_NB_CA_diangle = -116.094
        self.CD1_CG_NB_CA_diangle1 = -121.202

        self.CG_NB_length = 1.54
        self.CG_NB_CA_angle = 117.651
        self.CG_NB_CA_C_diangle = 88.315
        self.CG_NB_CA_C_diangle1 = 86.201

        self.a = -121.093
        self.a1 = -116.058
        self.b = 59.0
        self.c = -137.9
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.536
        self.C_N_CA_angle = 35.003

        self.CA_C_N_CD1_diangle = -173.574
        self.CA_C_N_CD1_diangle1 = -171.269

        # secondleft
        self.N_CL_length = 1.51
        self.N_CL_C15_angle = 106.96
        self.N_CL_C15_C14_diangle = -51.84

        self.CD1_N_CL_angle = 110.66
        self.CD1_N_CL_C15_diangle = -54.47

        self.CG_CD1_N_CL_diangle = -54.4

        #sidechain
        self.S_NB_length = 1.90
        self.S_NB_CG_angle = 109.66
        self.S_NB_CG_CD1_diangle = 103.41

        self.O1_S_length = 1.90
        self.O1_S_NB_angle = 109.06
        self.O1_S_NB_CG_diangle = -56.62

        self.O2_S_length = 1.90
        self.O2_S_NB_angle = 109.33
        self.O2_S_NB_CA_diangle = 55.75

        self.residue_name = "L"

    def inputRotamers(self, rotamers: List[float]) -> None:
        try:
            self.N_CA_CB_CG_diangle = rotamers[0]
            self.CA_CB_CG_CD1_diangle = rotamers[1]
            if self.CA_CB_CG_CD1_diangle > 0:
                self.CA_CB_CG_CD2_diangle = rotamers[1] - 180.0
            else:
                self.CA_CB_CG_CD2_diangle = rotamers[1] + 180.0
        except IndexError:
            print("Input Rotamers List: not long enough")
            self.N_CA_CB_CG_diangle = -66.4
            self.CA_CB_CG_CD1_diangle = 96.3
            self.CA_CB_CG_CD2_diangle = self.CA_CB_CG_CD1_diangle - 180.0

#helix_second_sidechain
class HssGeo(Geo):
    def __init__(self):
        self.CA_NB_length = 1.46
        self.CA_C_length = 1.52
        self.NB_CA_C_angle = 111.425

        self.CA_C_N_angle = 114.536
        self.NB_CA_C_N_diangle = -146.253
        self.NB_CA_C_N_diangle1 = -142.754
        self.C_N_CD1_angle = 120.434

        self.C_O_length = 1.225
        self.CA_C_O_angle = 119.92
        self.NB_CA_C_O_diangle = 41.184
        self.NB_CA_C_O_diangle1 = 44.937

        self.N_CD1_length = 1.427
        self.N_CD1_CG_angle = 112.261
        self.N_CD1_CG_NB_diangle = 61.910
        self.N_CD1_CG_NB_diangle1 = 66.252

        self.CD1_CG_length = 1.45
        self.CD1_CG_NB_angle = 114.475
        self.CD1_CG_NB_CA_diangle = -116.094
        self.CD1_CG_NB_CA_diangle1 = -121.202

        self.CG_NB_length = 1.54
        self.CG_NB_CA_angle = 117.651
        self.CG_NB_CA_C_diangle = 88.315
        self.CG_NB_CA_C_diangle1 = 86.201

        self.a = -121.093
        self.a1 = -116.058
        self.b = 59.0
        self.c = -137.9
        self.peptide_bond = 1.33
        self.CA_C_N_angle = 114.536
        self.C_N_CA_angle = 35.003

        self.CA_C_N_CD1_diangle = -173.574
        self.CA_C_N_CD1_diangle1 = -171.269

        # secondleft
        self.N_CL_length = 1.51
        self.N_CL_C15_angle = 106.96
        self.N_CL_C15_C14_diangle = -51.84

        self.CD1_N_CL_angle = 110.66
        self.CD1_N_CL_C15_diangle = -54.47

        self.CG_CD1_N_CL_diangle = -54.4

        # sidechain
        self.S_NB_length = 1.90
        self.S_NB_CG_angle = 108.29
        self.S_NB_CG_CD1_diangle = 86.41

        self.O1_S_length = 1.90
        self.O1_S_NB_angle = 107.64
        self.O1_S_NB_CA_diangle = 51.51

        self.O2_S_length = 1.90
        self.O2_S_NB_angle = 106.09
        self.O2_S_NB_CG_diangle = -42.32

        self.C1_S_length = 1.81
        self.C1_S_NB_angle = 108.85
        self.C1_S_NB_CG_diangle = 74.19

        self.CZ_C1_length = 1.54
        self.CZ_C1_S_angle = 109.54
        self.CZ_C1_S_NB_diangle = -135.16

        self.NH_CZ_length = 1.51
        self.NH_CZ_C1_angle = 109.38
        self.NH_CZ_C1_S_diangle = -179.96

        self.residue_name = "B"

#loop_first_sidechain
class LfsGeo(Geo):
    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98

        # sidechain
        self.C9_N2_length = 1.517
        self.C9_N2_C4_angle = 101.04
        self.C9_N2_C4_C3_diangle = 177.02

        self.C10_C9_length = 1.54
        self.C10_C9_N2_angle = 105.60
        self.C10_C9_N2_C4_diangle = -108.16

        self.C11_C10_length = 1.535
        self.C11_C10_C9_angle = 105.60
        self.C11_C10_C9_N2_diangle = 29.75

        self.residue_name = "C"

#loop_second_sidechain
class LssGeo(Geo):
    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98

        # sidechain

        self.residue_name = "D"

#loop_third_sidechain
class LtsGeo(Geo):
    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98

        # sidechain

        self.residue_name = "E"

#loop_fourth_sidechain
class LfoGeo(Geo):
    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98

        # sidechain

        self.residue_name = "F"

#loop_fifth_sidechain
class LfiGeo(Geo):
    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98

        # sidechain

        self.residue_name = "G"
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
        self.NB_CA_C_O_diangle=41.184
        self.NB_CA_C_O_diangle1 =44.937

        self.N_CD1_length=1.427
        self.N_CD1_CG_angle=112.261
        self.N_CD1_CG_NB_diangle=66.18
        self.N_CD1_CG_NB_diangle1=66.252

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

class Linker1Geo(Geo):

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

class Linker2Geo(Geo):

    def __init__(self):
        self.NL_C_length = 1.51
        self.NL_C_CA_angle = 105.23
        self.NL_C_CA_NB_diangle = -178.01

        self.C1_NL_length = 1.51
        self.C1_NL_C_angle = 113.89
        self.C1_NL_C_CA_diangle = -177.77

        self.C2_C1_length = 1.54
        self.C2_C1_NL_angle = 104.31
        self.C2_C1_NL_C_diangle = 178.46

        self.C3_C2_length = 1.54
        self.C3_C2_C1_angle = 106.19
        self.C3_C2_C1_NL_diangle = -57.36

        self.C4_C3_length = 1.54
        self.C4_C3_C2_angle = 111.33
        self.C4_C3_C2_C1_diangle = 170.82

        self.N2_C4_length = 1.51
        self.N2_C4_C3_angle = 109.10
        self.N2_C4_C3_C2_diangle = -175.94

        self.C5_N2_length = 1.51
        self.C5_N2_C4_angle = 105.06
        self.C5_N2_C4_C3_diangle = 53.18

        self.C6_C5_length = 1.54
        self.C6_C5_N2_angle = 109.58
        self.C6_C5_N2_C4_diangle = -173.81

        self.N3_C6_length = 1.50
        self.N3_C6_C5_angle = 111.23
        self.N3_C6_C5_N2_diangle = 53.68

        self.C7_N3_length = 1.51
        self.C7_N3_C6_angle = 103.51
        self.C7_N3_C6_C5_diangle = -163.88

        self.C8_C7_length = 1.54
        self.C8_C7_N3_angle = 112.26
        self.C8_C7_N3_C6_diangle = 164.30

        self.N4_C8_length = 1.50
        self.N4_C8_C7_angle = 108.25
        self.N4_C8_C7_N3_diangle = 70.50

        self.C13_N4_length = 1.505
        self.C13_N4_C8_angle = 110.62
        self.C13_N4_C8_C7_diangle = -164.22

        self.C14_C13_length = 1.534
        self.C14_C13_N4_angle = 104.67
        self.C14_C13_N4_C8_diangle = 176.52

        self.C15_C14_length = 1.533
        self.C15_C14_C13_angle = 111.43
        self.C15_C14_C13_N4_diangle = -171.23

        self.CL_C15_length = 1.534
        self.CL_C15_C14_angle = 103.24
        self.CL_C15_C14_C13_diangle = 175.98


        self.residue_name = "J"


def geometry(AA: str) -> Geo:
    """Generates the geometry of the requested amino acid.
    The amino acid needs to be specified by its single-letter
    code. If an invalid code is specified, the function
    returns the geometry of Glycine."""

    if AA == "L":
        return HfsGeo()
    elif AA == "B":
        return HssGeo()
    elif AA == "C":
        return LfsGeo()
    elif AA == "D":
        return LssGeo()
    elif AA == "E":
        return LtsGeo()
    elif AA == "F":
        return LfoGeo()
    elif AA == "G":
        return LfiGeo()
    elif AA == "H":
        return AAGeo()
    elif AA == "I":
        return Linker1Geo()
    elif AA == "A":
        return AlaGeo()
    else:
        return Linker2Geo()
