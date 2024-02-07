import lhapdf
import scipy
import numpy as np

#Pre-calculation variables:
s = 1000
M = 800
g = 1
scale = s


DCASIMIR = 1
SUMOVERREPS = 1

#Fermionic to start with.

#Import the gluon PDF set. Decide on the most appropriate one later.
p = lhapdf.mkPDF("CT10nlo", 0)
p = lhapdf.mkPDF("CT10nlo/0")

#Need to test the PDF set does as expected.
f_G = lambda x: p.xfxQ2(21, x, s)

##Partonic cross section. 
#Gluon fusion matrix element over s (i.e. eq 70 in draft)
MatrixElement_s = - Ga * Gb * SUMOVERREPS
#Partonic cross section
sigma_hat = - (MatrixElement_s / s) * 1/(12 * np.pi) * (1 + 2 M**2 / s) * np.sqrt(1 - 4*M**2 / s) * DCASIMIR * g**4 / 4


#Computing the convolution of the gluon PDFs with the partonic cross section (i.e. eq 71 & 72 in draft)
partonic = lambda x, y: f_G(x) * f_G(y) * sigma_hat(x*y*s) / 8
hadronic, err = scipy.integrate.quad(partonic, 0, np.infty)

print(hadronic)
print(err)