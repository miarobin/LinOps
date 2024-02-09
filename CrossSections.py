import lhapdf
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

'''
TO DO:
- Sort out alpha_s.
- Constants e.g. M, g, ...
- Add in other processes (e.g. non-fermions, different production mechanisms, ...).
- Check sigma_hat (i.e. manually set range of validity).
- Check against a known process.
'''

#Pre-calculation variables:
M = 250 #GeV
g = 1
QUARKS = {'u':[+2/3,lambda x: lhapdf.mkPDF("CT10nlo", 2)], 
          'd':[-1/3,lambda x: lhapdf.mkPDF("CT10nlo", 1)], 
          's':[-1/3,lambda x: lhapdf.mkPDF("CT10nlo", 3)]}
ANTIQUARKS = {'u':[-2/3,lambda x: lhapdf.mkPDF("CT10nlo", -2)], 
              'd':[+1/3,lambda x: lhapdf.mkPDF("CT10nlo", -1)], 
              's':[+1/3,lambda x: lhapdf.mkPDF("CT10nlo", -3)]}

Ga = Gb = 1
DCASIMIR = 1
SUMOVERREPS = 1

#Fermionic to start with.

#Import the gluon PDF set. Decide on the most appropriate one later.
p = lhapdf.mkPDF("CT10nlo", 0)
p = lhapdf.mkPDF("CT10nlo/0")

#This probably needs mofidying for new particles.
alpha_S = lhapdf.mkAlphaS("CT10nlo")


#Need to test the PDF set does as expected.
f_G = lambda x, s_max: p.xfxQ2(21, x, s_max)


##Test of SM Drell-Yan.
test_s=np.power(np.arange(0,1000),2) #GeV
result = []
for ts in test_s:
    sigma_partonic_DY = lambda x1, x2: (4*np.pi*alpha_S.alphasQ(ts)**2 / (ts**2))* np.sum([x1*QUARKS[q][1](x1)*x2*ANTIQUARKS[q][1](x2) * QUARKS[q][0]**2 for q in ['u','d','s']]) 
    sigma_hadronic_DY, err = integrate.nquad(sigma_partonic_DY, [[0, 1],[0, 1]])
    result.append(sigma_hadronic_DY)
plt.plot(np.power(test_s,0.5), result)
plt.xlabel("Centre of Mass Energy")
plt.ylabel("Cross Section Drell Yan")
plt.savefig("DrellYanTest.pdf", format="pdf", bbox_inches="tight")


##Partonic cross section. 
#Gluon fusion matrix element over s (i.e. eq 70 in draft)
MatrixElement_s = - SUMOVERREPS
#Partonic cross section
F_fermion = lambda s: 2 * np.sqrt(1-4*M**2/s) * (1+2*M**2/s) if s>=4*M**2 else 0
F_scalar = lambda s: (1-4*M**2)**(3/2) if s>=4*M**2 else 0


#Run over s_max:
s_maxs = np.power(np.arange(100,3000),2)
hadronics = []
for s_max in s_maxs:
    #Computing the convolution of the gluon PDFs with the partonic cross section (i.e. eq 71 & 72 in draft)(NOT FINISHED)
    
    #Constants in front of the integral:
    consts = np.pi * alpha_S.alphasQ(s_max)**2
    #Perform the convolution integral:
    partonic = lambda x, y: (f_G(x,s_max)*f_G(y,s_max)*F_fermion(x*y*s_max) / (x*y)) * consts
    hadronic, err = integrate.nquad(partonic, [[0, 1],[0, 1]])/s_max
    hadronics.append(hadronic)


plt.plot(np.power(s_maxs,0.5), hadronics)
plt.xlabel("Centre of Mass Energy")
plt.ylabel("Cross Section")
plt.savefig("testing.pdf", format="pdf", bbox_inches="tight")

    