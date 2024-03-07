import lhapdf
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt


#Import the PDF sets. Decide on the most appropriate one later.
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
#Strong coupling constant running.
alpha_S = lhapdf.mkAlphaS("CT18NNLO")


plt.figure()
#NEW STUFF
#Partonic cross sections:
def G_scalar(betasq):
    if betasq>=0:
        beta = np.sqrt(betasq)
        r=1
        return beta*(2-betasq)/2 - (1-betasq**2)*np.arctanh(beta)/2 +\
            r*(beta*(3-5*betasq)/24 - (1-betasq)**2*np.arctanh(beta)/4)
    else:
        return 0

def G_fermion(betasq):
    if betasq>=0:
        beta = np.sqrt(betasq)
        r=1
        return beta*(betasq-2) + (3-betasq**2)*np.arctanh(beta)+\
            r*(5*beta*(betasq-3)/12 + (1-betasq)**2*np.arctanh(beta)/2)
    else:
        return 0
    

def F_scalar(betasq):
    if betasq>=0:
        beta = np.sqrt(betasq)

        return beta**3/3
    else:
        return 0
    
def F_fermion(betasq):
    if betasq>=0:
        beta = np.sqrt(betasq)
        return (2/3)*beta*(3-betasq)
    else:
        return 0
        

'''
#HADRONIC CROSS-SECTION
Mnews=np.linspace(500,1000,num=3)
LHC = (13.5e3)**2 #GeV^2

print(alpha_S.alphasQ2(100**2))
print(alpha_S.alphasQ2(500**2))

results = []
for Mn in Mnews:
    betasq = lambda x,y: 1 - 4*Mn**2/(LHC*x*y)
    constsGG=1


    sigma_GGf = integrate.nquad(lambda x,y: f_G(x,LHC)*f_G(y,LHC)*G_fermion(betasq(x,y))/(x*y)**2,[[0.001,1],[0.001,1]])[0]
    sigma_GGs = integrate.nquad(lambda x,y: f_G(x,LHC)*f_G(y,LHC)*G_scalar(betasq(x,y))/(x*y)**2,[[0.001,1],[0.001,1]])[0]

    
    #Notice we've calculated Q_Y,q (hypercharges) here & summed over left & right charges.
    sigma_qqYs = integrate.nquad(lambda x,y: F_scalar(betasq(x,y))/(x*y)**2 *\
                                np.sum([(QUARKS[q][1](x,LHC)*ANTIQUARKS[q][1](y,LHC))*(QUARKS[q][0]**2+(1/6)**2) for q in ['u','d']]),[[0.001,1],[0.001,1]])[0]
    
    sigma_qqYf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                                np.sum([(QUARKS[q][1](x,LHC)*ANTIQUARKS[q][1](y,LHC))*(QUARKS[q][0]**2+(1/6)**2) for q in ['u','d']]),[[0.001,1],[0.001,1]])[0]
    

    sigma_qqLf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                                (2*(QUARKS['u'][1](x,LHC)*ANTIQUARKS['d'][1](y,LHC) + QUARKS['d'][1](x,LHC)*ANTIQUARKS['u'][1](y,LHC)) +\
                                (QUARKS['u'][1](x,LHC)*ANTIQUARKS['u'][1](y,LHC) + QUARKS['d'][1](x,LHC)*ANTIQUARKS['d'][1](y,LHC))),[[0.001,1],[0.001,1]])[0]/4


    sigma_qqLs = integrate.nquad(lambda x,y: F_scalar(betasq(x,y))/(x*y)**2 *\
                                (2*(QUARKS['u'][1](x,LHC)*ANTIQUARKS['d'][1](y,LHC) + QUARKS['d'][1](x,LHC)*ANTIQUARKS['u'][1](y,LHC)) +\
                                (QUARKS['u'][1](x,LHC)*ANTIQUARKS['u'][1](y,LHC) + QUARKS['d'][1](x,LHC)*ANTIQUARKS['d'][1](y,LHC))),[[0.001,1],[0.001,1]])[0]/4

    results.append([sigma_GGf,sigma_GGs,sigma_qqYs,sigma_qqYf,sigma_qqLf,sigma_qqLs])

results = np.array(results)
plt.plot(Mnews,results[:,0],color='red')
plt.plot(Mnews,results[:,1],color='blue')
plt.plot(Mnews,results[:,2],color='black')
plt.plot(Mnews,results[:,3],color='grey')
plt.plot(Mnews,results[:,4],color='orange')
plt.plot(Mnews,results[:,5],color='green')
plt.xlabel("Mass of New Particle")
plt.ylabel("Cross Section")
plt.savefig("testing.pdf", format="pdf", bbox_inches="tight")'''

plt.figure()
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
def dG_dt(s,M,t):
    R = 1
    return (2/s) * ( (1/4)*(s**2/(t*(-t-s)) - 2) + (M**2*s/((-s-t)*t)) * (1 - M**2*s/((-t-s)*t)) \
                                        + R*( (1/2)*(((-t-s)*t)/s**2 - 1) + (M**2/s)*(s*M**2/((-s-t)*t) - 1)))
def dGs_dx(x,beta):
    R = 1
    return (1/2) + ((beta**2-1)/4)* 1/(x**2 + x) * (((beta**2-1)/4)* 1/(x**2 + x) + 1) + R*( (1/8) + (x**2 + x)/2 + ((beta**2-1)/4)*(((beta**2-1)/4)* 1/(x**2 + x) + 1))
    
#Use beta = 1/2 and M=1000GeV.
xs = np.linspace(-(1+1/2)/2, -(1-1/2)/2)
plt.plot(xs,4*1000**2/(1-.25)*dG_dt(4*1000**2/(1-.25),1000,xs*4*1000**2/(1-.25)))
plt.savefig("dGdtTest.pdf", format="pdf", bbox_inches="tight")

plt.figure()
#Use beta = 1/2 and M=1000GeV.
xs = np.linspace(-(1+1/2)/2, -(1-1/2)/2)
plt.plot(xs,dGs_dx(xs,1/2))
plt.savefig("dGsdxTest.pdf", format="pdf", bbox_inches="tight")


#3: Test against Rodrigo's partonic cross sections (gluonic and fermionic)

resultsG = []; resultsF = []; betasqs = np.linspace(0,1)
for betasq in betasqs:
    resultsG.append([G_fermion(betasq),G_scalar(betasq)])
    resultsF.append([F_fermion(betasq),F_scalar(betasq)])
resultsG = np.array(resultsG); resultsF = np.array(resultsF)
xs=((1-betasqs)/4)**0.5

plt.figure()
plt.plot(xs,resultsG[:,0],color='red')
plt.plot(xs,resultsG[:,1],color='blue')

plt.savefig("partonicG.pdf", format="pdf", bbox_inches="tight")

plt.figure()
plt.plot(xs,resultsF[:,0]/2,color='red')
plt.plot(xs,resultsF[:,1]/2,color='blue')

plt.savefig("partonicF.pdf", format="pdf", bbox_inches="tight")