import csv
from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import scipy.interpolate

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["font.size"]= 13


#Read into a numpy array the edited data.
delimiter = ','
dataF = np.array(np.genfromtxt('HadronicCrossSecFermions.csv', delimiter=delimiter, skip_header=1, dtype=None))
dataS = np.array(np.genfromtxt('HadronicCrossSecScalar.csv', delimiter=delimiter, skip_header=1, dtype=None))

Mnews = np.linspace(500,2500,num=10)

NPf = []; NPs = []
for item in dataS:
	NPs.append(item)
for item in dataF:
	NPf.append(item)
NPf = np.array(NPf); NPs = np.array(NPs)

fig, ax = plt.subplots()

#Bounds on Plot



#Interpolate so that we can draw the bound better
MnewsPrime = np.linspace(500,2500,num=100)
interpolatedNPf = scipy.interpolate.interp1d(Mnews, np.log(NPf[:,0]), kind='linear')

allowedThetaFe = [i for i, m in enumerate(MnewsPrime) if not 50<=m<=600]
disallowedThetaFe = [i for i, m in enumerate(MnewsPrime) if 50<=m<=600] + [allowedThetaFe[0]]

ax.plot(MnewsPrime[allowedThetaFe],np.exp(interpolatedNPf(MnewsPrime[allowedThetaFe])),color='red'); ax.plot(Mnews,NPs[:,0],color='red',linestyle='dashed')
ax.plot(MnewsPrime[disallowedThetaFe],np.exp(interpolatedNPf(MnewsPrime[disallowedThetaFe])),color='red',linestyle='dotted')


ax.plot(Mnews,NPf[:,1],color='blue'); ax.plot(Mnews,NPs[:,1],color='blue',linestyle='dashed')
ax.plot(Mnews,NPf[:,2],color='black'); ax.plot(Mnews,NPs[:,2],color='black',linestyle='dashed')
ax.plot(Mnews,NPf[:,3],color='grey'); ax.plot(Mnews,NPs[:,3],color='grey',linestyle='dashed')
ax.plot(Mnews,NPf[:,4],color='orange'); ax.plot(Mnews,NPs[:,4],color='orange',linestyle='dashed')
ax.plot(Mnews,NPf[:,5],color='green'); ax.plot(Mnews,NPs[:,5],color='green',linestyle='dashed')
ax.plot(Mnews,NPf[:,4]-NPf[:,1],color='purple'); ax.plot(Mnews,NPs[:,4]-NPs[:,1],color='purple',linestyle='dashed')
ax.set_xscale('log')
ax.set_yscale('log')


ax.get_xaxis().set_major_formatter(matplotlib.ticker.LogFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.LogFormatter())
ax.set_xticks([500,1000,2000])

ax.set_xlabel("Mass of New Particle, M (GeV)")
ax.set_ylabel("Cross Section (pb)")

#plt.show()
plt.savefig("CrossSecPlot_CMSBounds.pdf", format="pdf", bbox_inches="tight")