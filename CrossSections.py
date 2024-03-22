import lhapdf
import numpy as np
from scipy import integrate
import csv
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.size"]= 12


def save_arrays_to_csv(file_path, column_titles, *arrays):
    # Transpose the arrays to align them by columns
    transposed_arrays = np.array(list(zip(*arrays)))
    

    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        #Writes the column titles.
        writer.writerow(column_titles)
        
        #Writes in data rows.
        for row in transposed_arrays:
            writer.writerow(row)

#Import the PDF sets. Decide on the most appropriate one later.
p = lhapdf.mkPDF("PDF4LHC15_nlo_100", 0)
p = lhapdf.mkPDF("PDF4LHC15_nlo_100/0")

#SM Particles [0]: LH Hypercharge, [1]: RH Hypercharge, [2]: Proton PDF
QUARKS = {'u':[+1/6,+2/3,lambda x,s_max: p.xfxQ2(2, x, s_max)], 
          'd':[+1/6,-1/3,lambda x,s_max: p.xfxQ2(1, x, s_max)]}
ANTIQUARKS = {'u':[-1/6,-2/3,lambda x,s_max: p.xfxQ2(-2, x, s_max)], 
              'd':[-1/6,+1/3,lambda x,s_max: p.xfxQ2(-1, x, s_max)]}

#New Particles [0]: Fundamental SU(2) rep., [1]: Fundamental colour rep., [2]: Hypercharge
Zs = {'Xi':[0, 0, 1/2],
      'Lambda': [1, 0, 0],
      'Omega': [0, 1, 1/6],
      'Sigma': [0, 0, 1/3],
      'Delta': [1, 0, 1/6],
      'Theta': [0, 1, 0]}


#Gluon PDF
f_G = lambda x, s_max: p.xfxQ2(21, x, s_max)
#Strong coupling constant running.
alpha_S = lhapdf.mkAlphaS("PDF4LHC15_nlo_100")
mW = 80; mZ = 91; v = 246
alpha_Y = ((2*mZ/v)**2 - (2*mW/v)**2)/(4*np.pi)
alpha_w = (2*mW/v)**2/(4*np.pi)


plt.figure()
#NEW STUFFO
#Partonic cross sections:
def G_scalar(betasq,r):
    if betasq>=0:
        beta = np.sqrt(betasq)
        return beta*(2-betasq)/2 - (1-betasq**2)*np.arctanh(beta)/2 +\
            r*(beta*(3-5*betasq)/24 - (1-betasq)**2*np.arctanh(beta)/4)
    else:
        return 0

def G_fermion(betasq,r):
    if betasq>=0:
        beta = np.sqrt(betasq)
        return beta*(betasq-2) + (3-betasq**2)*np.arctanh(beta)+\
            r*(beta*(5*betasq-9)/12 + (1-betasq)**2*np.arctanh(beta)/2)
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
        


#HADRONIC CROSS-SECTION
Mnews=np.linspace(500,2500,num=10)
LHC = (13.5e3)**2 #GeV^2

print(f"Alpha_s at Z pole: {alpha_S.alphasQ2(90**2)}")
print(f"Alpha_s at 500 GeV: {alpha_S.alphasQ2(500**2)}")

NPs = []; NPf = []
for Mn in Mnews:
    betasq = lambda x,y: 1 - 4*Mn**2/(LHC*x*y)
    
    dL = lambda n: n+1 ; DL = lambda n: n*(n+1)*(n+2)/(3*2*2)
    dc = lambda n,m: (m+1)*(n+1)*(n+m+2)/2 ; Dc = lambda n,m :(m**3 + n**3 + 3*(n+m) + m*n) * dc(n,m)/ (4*3*2)
    
    constsGG = lambda nL, nC: np.pi*alpha_S.alphasQ2(91**2)**2 * dL(nL) * 4 * Dc(nC,0)**2 / (LHC*dc(nC,0))
    constsqqL = lambda nL, nC: np.pi*alpha_w**2 * dc(nC,0) * DL(nL) / (6*LHC)
    constsqqY = lambda nL, nC, QY: np.pi * alpha_Y**2 * QY**2 * dc(nC,0) * dL(nL) / (6*LHC)

    #NOTE the PDFs from LHAPDF are of the form f_LHAPDF = xf_DRAFT(x).

    #QUARKS These are the same for any new particle so just do it once.
    sigma_qqYs = integrate.nquad(lambda x,y: F_scalar(betasq(x,y))/(x*y)**2 *\
                            np.sum([(QUARKS[q][2](x,x*y*LHC)*ANTIQUARKS[q][2](y,x*y*LHC))*(QUARKS[q][0]**2/2+QUARKS[q][1]**2/2) for q in ['u','d']]),[[0.001,1],[0.001,1]])[0]
    

    sigma_qqYf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                            np.sum([(QUARKS[q][2](x,x*y*LHC)*ANTIQUARKS[q][2](y,x*y*LHC))*(QUARKS[q][0]**2/2+QUARKS[q][1]**2/2) for q in ['u','d']]),[[0.001,1],[0.001,1]])[0]
    
    #Note for dividing by 16, 2x2 is from Pauli matrices/2 and 2x2 from only taking LH quarks.
    sigma_qqLf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                            (2*(QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC)) +\
                            (QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC))),[[0.001,1],[0.001,1]])[0]/16


    sigma_qqLs = integrate.nquad(lambda x,y: F_scalar(betasq(x,y))/(x*y)**2 *\
                            (2*(QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC)) +\
                            (QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC))),[[0.001,1],[0.001,1]])[0]/16
    
    #GLUONS These differ depending on the new particle so have to do multiple times.
    r = lambda n, m: 3 * dc(n,m) / (Dc(n,m) * (3**2-1))

    fresults_ = []; sresults_ = []
    for Z in Zs.keys():
        if Zs[Z][1]!= 0:
            sigma_GGf = integrate.nquad(lambda x,y: f_G(x,x*y*LHC)*f_G(y,x*y*LHC)*G_fermion(betasq(x,y),r(Zs[Z][1],0))/(x*y)**2,[[0.001,1],[0.001,1]])[0]
            sigma_GGs = integrate.nquad(lambda x,y: f_G(x,x*y*LHC)*f_G(y,x*y*LHC)*G_scalar(betasq(x,y),r(Zs[Z][1],0))/(x*y)**2,[[0.001,1],[0.001,1]])[0]
        else: sigma_GGf = 0; sigma_GGs = 0
        ZF = sigma_GGf*constsGG(*Zs[Z][0:2]) + sigma_qqYf*constsqqY(*Zs[Z]) + sigma_qqLf*constsqqL(*Zs[Z][0:2])
        ZS = sigma_GGs*constsGG(*Zs[Z][0:2]) + sigma_qqYs*constsqqY(*Zs[Z]) + sigma_qqLs*constsqqL(*Zs[Z][0:2])
        fresults_.append(ZF)
        sresults_.append(ZS)
    NPf.append(fresults_)
    NPs.append(sresults_)

NPf = np.array(NPf)/(2.56819e-9); NPs = np.array(NPs)/(2.56819e-9)
plt.plot(Mnews,NPf[:,0],color='red'); plt.plot(Mnews,NPs[:,0],color='red',linestyle='dashed')
plt.plot(Mnews,NPf[:,1],color='blue'); plt.plot(Mnews,NPs[:,1],color='blue',linestyle='dashed')
plt.plot(Mnews,NPf[:,2],color='black'); plt.plot(Mnews,NPs[:,2],color='black',linestyle='dashed')
plt.plot(Mnews,NPf[:,3],color='grey'); plt.plot(Mnews,NPs[:,3],color='grey',linestyle='dashed')
plt.plot(Mnews,NPf[:,4],color='orange'); plt.plot(Mnews,NPs[:,4],color='orange',linestyle='dashed')
plt.plot(Mnews,NPf[:,5],color='green'); plt.plot(Mnews,NPs[:,5],color='green',linestyle='dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Mass of New Particle (GeV)")
plt.ylabel("Cross Section (pb)")
plt.savefig("testing.pdf", format="pdf", bbox_inches="tight")

#FILE WRITER
#Column Titles
column_titles = ['Xi', 'Lambda', 'Omega', 'Simga', 'Delta', 'Theta']
# File path to save the CSV
file_path_f = f'HadronicCrossSecFermions.csv'
file_path_s = f'HadronicCrossSecScalar.csv'
save_arrays_to_csv(file_path_f, column_titles, NPf[:,0],NPf[:,1],NPf[:,2],NPf[:,3],NPf[:,4],NPf[:,5])
save_arrays_to_csv(file_path_s, column_titles, NPs[:,0],NPs[:,1],NPs[:,2],NPs[:,3],NPs[:,4],NPs[:,5])
plt.figure()

#TESTS

#1: pdf test using SM Drell-Yan
MDYs=np.linspace(100,400,num=50) #GeV
result = []
for MDY in MDYs:
    tau = MDY**2/(1800)**2
    def d_sigma_dMdY(x):
        return (8*np.pi*(7.297e-3)**2 / (3*3*MDY**3)) * \
                np.sum([(QUARKS[q][2](np.sqrt(tau)*np.exp(x),MDY**2)*QUARKS[q][2](np.sqrt(tau)*np.exp(-x),MDY**2) + 
                         ANTIQUARKS[q][2](np.sqrt(tau)*np.exp(x),MDY**2)*ANTIQUARKS[q][2](np.sqrt(tau)*np.exp(-x),MDY**2))* QUARKS[q][0]**2 for q in ['u','d']])

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
    resultsG.append([G_fermion(betasq,1),G_scalar(betasq,1)])
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


#4: Testing against Leptoquarks at the Tevatron
'''
MLQs=np.linspace(100,1400,num=5)
TEV = (1.8e3)**2 #GeV^2

results = []
for Mn in MLQs:
    betasq = lambda x,y: 1 - 4*Mn**2/(TEV*x*y)
    
    dL = lambda n: n+1 ; DL = lambda n: n*(n+1)*(n+2)/(3*2*2)
    dc = lambda n,m: (m+1)*(n+1)*(n+m+2)/2 ; Dc = lambda n,m :(m**3 + n**3 + 3*(n+m) + m*n) * dc(n,m)/ (4*3*2)
    
    constsGG = lambda nL, nC: np.pi * dL(nL) * 4 * Dc(nC,0)**2 / (TEV*dc(nC,0))
    constsqqL = lambda nL, nC: np.pi*alpha_w**2 * dc(nC,0) * dL(nL) / TEV
    constsqqY = lambda nL, nC, QY: np.pi * alpha_Y**2 * QY**2 * dc(nC,0) * dL(nL) / TEV

    #NOTE the PDFs from LHAPDF are of the form f_LHAPDF = xf_DRAFT(x).
    #GLUONS.
    r = lambda n, m: 3 * dc(n,m) / (Dc(n,m) * (3**2-1))
    sigma_GGs = integrate.nquad(lambda x,y: alpha_S.alphasQ2(91**2)**2*f_G(x,(x*y*TEV))*f_G(y,(x*y*TEV))*G_scalar(betasq(x,y),r(Zs['Theta'][1],0))/(x*y)**2,[[0.001,1],[0.001,1]])[0]

    print(r(Zs['Theta'][1],0))
    #Multiply by constants & add to result array.
    LQGG = sigma_GGs*constsGG(*Zs['Theta'][0:2]) 
    
    results.append(LQGG)

results = np.array(results)
plt.figure()

plt.plot(MLQs,results/(2.56819e-9),color='red')
plt.yscale('log')
plt.xlabel("Mass of New Particle")
plt.ylabel("Cross Section")
plt.savefig("LeptoquarkTest.pdf", format="pdf", bbox_inches="tight")

plt.figure()
'''

#Testing against heavy leptons at LHC
MLQs=np.linspace(500,2000,num=5)
LHC = (14e3)**2 #GeV^2

results = []
for Mn in MLQs:
    betasq = lambda x,y: 1 - 4*Mn**2/(LHC*x*y)
    
    dL = lambda n: n+1 ; DL = lambda n: n*(n+1)*(n+2)/(3*2*2)
    dc = lambda n,m: (m+1)*(n+1)*(n+m+2)/2 ; Dc = lambda n,m :(m**3 + n**3 + 3*(n+m) + m*n) * dc(n,m)/ (4*3*2)
    
    constsqqL = lambda nL, nC: np.pi * alpha_w**2 * dc(nC,0) * DL(nL) / (6*LHC)
    constsqqY = lambda nL, nC, QY: np.pi * alpha_Y**2 * QY**2 * dc(nC,0) * dL(nL) / (6*LHC)

    #NOTE the PDFs from LHAPDF are of the form f_LHAPDF = xf_DRAFT(x).

    sigma_qqYf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                            np.sum([(QUARKS[q][2](x,x*y*LHC)*ANTIQUARKS[q][2](y,x*y*LHC))*(QUARKS[q][0]**2/2+QUARKS[q][1]**2/2) for q in ['u','d']]),[[0.001,1],[0.001,1]])[0]
    
    sigma_qqLf = integrate.nquad(lambda x,y: F_fermion(betasq(x,y))/(x*y)**2 *\
                            (2*(QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC)) +\
                            (QUARKS['u'][2](x,x*y*LHC)*ANTIQUARKS['u'][2](y,x*y*LHC) + QUARKS['d'][2](x,x*y*LHC)*ANTIQUARKS['d'][2](y,x*y*LHC))),[[0.001,1],[0.001,1]])[0]/16

    #Multiply by constants & add to result array.
    
    QQLL = constsqqL(1,0)*sigma_qqLf + constsqqY(1,0,-1)*sigma_qqYf
    
    results.append(QQLL)

results = np.array(results)
plt.figure()

plt.plot(MLQs,results/(2.56819e-12),color='red')
plt.yscale('log')
plt.xlabel("Mass of New Particle")
plt.ylabel("Cross Section")
plt.savefig("HeavyLeptonTest.pdf", format="pdf", bbox_inches="tight")
