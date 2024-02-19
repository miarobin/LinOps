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


plt.figure()



#NEW STUFF
#Partonic cross sections
def dG_dt(s,M,t):
    R = 1
    return (2/s) * ( (1/4)*(s**2/(t*(-t-s)) - 2) + (M**2*s/((-s-t)*t)) * (1 - M**2*s/((-t-s)*t)) \
                                        + R*( (1/2)*(((-t-s)*t)/s**2 - 1) + (M**2/s)*(s*M**2/((-s-t)*t) - 1)))

def G_fermion(s,M):
    betasq = 1 - 4*M**2/s
    if betasq >= 0:
        beta = np.sqrt(betasq)
        return integrate.quad(lambda t: dG_dt(s,M,t), -s*(1+beta)/2, -s*(1-beta)/2)[0]
    else:
        return 0
    
def dGs_dx(x,beta):
    R = 1
    return (1/2) + ((beta**2-1)/4)* 1/(x**2 + x) * (((beta**2-1)/4)* 1/(x**2 + x) + 1) + R*( (1/8) + (x**2 + x)/2 + ((beta**2-1)/4)*(((beta**2-1)/4)* 1/(x**2 + x) + 1))
    
def G_scalar(s,M):
    betasq = 1 - 4*M**2/s
    if betasq>=0:
        beta = np.sqrt(betasq)
        return integrate.quad(lambda x: dGs_dx(x,beta), -(1+beta)/2, -(1-beta)/2)
    else:
        return 0
    
def F_scalar(s,M):
    betasq = 1 - 4*M**2/s
    if betasq>=0:
        beta = np.sqrt(betasq)
        return beta**3/3
    else:
        return 0
    
def F_fermion(s,M):
    betasq = 1 - 4*M**2/s
    if betasq>=0:
        beta = np.sqrt(betasq)
        return (2/3)*beta*(1-beta**2)
        
    
#COLOUR CROSS-SECTION
Mnews=np.linspace(600,1000,num=50)
LHC = (13.5e3)**2 #GeV^2

result = []
for Mn in Mnews:
    consts=1
    #sigma_color = consts_color*integrate.nquad(lambda x,y: f_G(x,LHC)*f_G(y,LHC)*G_fermion(x*y*LHC,Mn)/(x*y)**2,[[0,1],[0,1]])
    sigma_qqY = integrate.nquad(lambda x,y: F_scalar(x*y*LHC,Mn)/(x*y)**2 *\
                                np.sum([QUARKS[q][1](x)*ANTIQUARKS[q][1](y)*QUARKS[q][0]**2 for q in ['u','d']]),[[0,1],[0,1]])
    result.append(sigma_qqY)

plt.plot(Mnews,result)
plt.xlabel("Mass of New Particle")
plt.ylabel("Cross Section")
plt.savefig("testing.pdf", format="pdf", bbox_inches="tight")


#TESTS

#1: pdf test using SM Drell-Yan
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

#2: Test of gluon partonic cross section against Rodrigo's
plt.figure()
#Use beta = 1/2 and M=1000GeV.
xs = np.linspace(-(1+1/2)/2, -(1-1/2)/2)
plt.plot(xs,4*1000**2/(1-.25)*dG_dt(4*1000**2/(1-.25),1000,xs*4*1000**2/(1-.25)))
plt.savefig("dGdtTest.pdf", format="pdf", bbox_inches="tight")

plt.figure()
#Use beta = 1/2 and M=1000GeV.
xs = np.linspace(-(1+1/2)/2, -(1-1/2)/2)
plt.plot(xs,dGs_dx(xs,1/2))
plt.savefig("dGsdxTest.pdf", format="pdf", bbox_inches="tight")