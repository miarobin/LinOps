import lhapdf
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt

'''
TO DO:
- Constants e.g. M, g, ...
- Add in other processes (e.g. non-fermions, different production mechanisms, ...).
'''


#Import the PDF set. Decide on the most appropriate one later.
p = lhapdf.mkPDF("CT18NNLO", 0)
p = lhapdf.mkPDF("CT18NNLO/0")
QUARKS = {'u':[+2/3,lambda x,s_max: p.xfxQ2(2, x, s_max)], 
          'd':[-1/3,lambda x,s_max: p.xfxQ2(1, x, s_max)], 
          's':[-1/3,lambda x,s_max: p.xfxQ2(3, x, s_max)],
          'c':[+2/3,lambda x,s_max: p.xfxQ2(4, x, s_max)]}
ANTIQUARKS = {'u':[-2/3,lambda x,s_max: p.xfxQ2(-2, x, s_max)], 
              'd':[+1/3,lambda x,s_max: p.xfxQ2(-1, x, s_max)], 
              's':[+1/3,lambda x,s_max: p.xfxQ2(-3, x, s_max)],
              'c':[-2/3,lambda x,s_max: p.xfxQ2(-4, x, s_max)]}

f_G = lambda x, s_max: p.xfxQ2(21, x, s_max)

alpha_S = lhapdf.mkAlphaS("CT18NNLO")



##Test of SM Drell-Yan.
MDYs=np.linspace(100,400,num=50) #GeV
result = []
for MDY in MDYs:
    tau = MDY**2/(1800)**2
    def d_sigma_dMdY(x):
        return (8*np.pi*(7.297e-3)**2 / (3*3*MDY**3)) * \
                np.sum([(QUARKS[q][1](np.sqrt(tau)*np.exp(x),MDY**2)*QUARKS[q][1](np.sqrt(tau)*np.exp(-x),MDY**2) + 
                         ANTIQUARKS[q][1](np.sqrt(tau)*np.exp(x),MDY**2)*ANTIQUARKS[q][1](np.sqrt(tau)*np.exp(-x),MDY**2))* QUARKS[q][0]**2 for q in ['u','d','s','c']])

    sigma_hadronic_DY, err = integrate.quad(d_sigma_dMdY, -1,1)
    result.append(sigma_hadronic_DY*0.5)

plt.plot(MDYs, np.array(result)*(0.3894*1e9))
plt.yscale('log')
plt.xlabel("Lepton invariant mass M (GeV)")
plt.ylabel(f"$d^2\sigma/dMdy$ for $|y|<1$ (pb/GeV)")
plt.savefig("DrellYanTest.pdf", format="pdf", bbox_inches="tight")

plt.figure()

#NEW STUFF
#Partonic cross sections
def G_gluon(s,M):
    R = 1
    if 1-4*M**2/s >= 0:
        integrand = lambda t1: (2/s) * ( -(1/4)*(s**2/(t1**2 + s*t1) + 2) - (M**2*s/(s**2*t1 + t1**2)) * (1 + M**2*s/(s*t1 + t1**2)) \
                                        + R*( -(1/2)*((s*t1 + t1**2)/s**2 + 1) - (M**2/s)*(s*M**2/(t1**2 + s*t1) + 1)))
        return integrate.quad(integrand, -s, +s)[0]
    else:
        return 0
    
#COLOUR CROSS-SECTION
Mnews=np.linspace(500,2000,num=50)
LHC = 100e3 #GeV

result = []
for Mn in Mnews:
    consts_color=1
    sigma_color = consts_color*integrate.nquad(lambda x,y: f_G(x,LHC)*f_G(y,LHC)*G_gluon(x*y*LHC,Mn)/(x*y)**2,[[0,1],[0,1]])

    result.append(sigma_color)


plt.plot(Mnews,result)
plt.xlabel("Mass of New Particle")
plt.ylabel("Cross Section")
plt.savefig("testing.pdf", format="pdf", bbox_inches="tight")
