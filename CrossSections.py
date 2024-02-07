import lhapdf
import numpy as np
from scipy import integrate

#Pre-calculation variables:
s_max = 1e4**2 #GeV
M = 100 #GeV
g = 1
scale = s_max

Ga = Gb = 1
DCASIMIR = 1
SUMOVERREPS = 1

#Fermionic to start with.

#Import the gluon PDF set. Decide on the most appropriate one later.
p = lhapdf.mkPDF("CT10nlo", 0)
p = lhapdf.mkPDF("CT10nlo/0")

alpha_S = lhapdf.mkAlphaS("CT10nlo")


#Need to test the PDF set does as expected.
f_G = lambda x: p.xfxQ2(21, x, s_max)

##Partonic cross section. 
#Gluon fusion matrix element over s (i.e. eq 70 in draft)
MatrixElement_s = - Ga * Gb * SUMOVERREPS
#Partonic cross section
sigma_hat = lambda s: - (MatrixElement_s / s) * 1/(12 * np.pi) * (1 + 2 * M**2 / s) * np.sqrt(1 - 4*M**2 / s) * DCASIMIR * g**4 / 4


#Computing the convolution of the gluon PDFs with the partonic cross section (i.e. eq 71 & 72 in draft)(NOT FINISHED)
partonic = lambda x, y: (f_G(x) * f_G(y) * sigma_hat(x*y*s_max) / (x*y)) * np.pi * alpha_S.alphasQ(x*y*s_max)**2 * DCASIMIR / (8 * DCASIMIR**2 * 12) 
hadronic, err = integrate.nquad(partonic, [[0, 1],[0, 1]])

print(hadronic)
print(err)