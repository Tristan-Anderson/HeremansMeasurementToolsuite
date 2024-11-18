import pandas as pd
import numpy as np
import pylab as pl
import utilities, windower
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit

### General function to extract carrier densities from measurements
### Selects single region
def quantum_hall(
		fn, # Df containing data
		r,  # Df (record)
		n,  # Measurement number
		x='B Field (T)', 
		y='P124A (V)', 
		f=None):
	""" 
	    Function which handles Hall Data. Requests singular file, with paramaters: x,y, and selx.
	    User can elect x-selection capabilities on call()
	"""
	e = 1.60217663*10**-19
	h = 6.62607015*10**-34
	Rk = h/(e**2)
	RName = "Klitzing Resistance"
	m = r[r["Measurement"]==n]
	I = m["I_ac"].values.tolist()[0]
	c = m["Correction"].values.tolist()[0]
	df = utilities.ReadDatFile(fn, correction=c, y=y)
	df[RName] = (df[y].values/I)/Rk
	dpts = []
	addition = [None,None]
	fig,ax = None,None
	while len(addition) == 2:
		addition = windower.selectWindower(df,x,\
				RName,clickpoints=True, fig=fig,\
				ax=ax,filename=fn)
		dpts += addition
	else:
		if len(addition) == 0:
			pass
		else:
			cols = [i for i in range(len(dpts))]
			results = {i:dpts[i] for i in cols}
			rdf = pd.DataFrame(results)
			print(rdf)
			with open("QH_Record_"+fn.split('.dat')[0]+'.csv','w') as recordfile:
				rdf.to_csv(recordfile)		
			plt.close()
		return dpts

def quantum_hall_density(
		fn,
		r,
		n,
		x='B Field (T)', 
		y='P124A (V)', 
		f=None):
	# elect peaks
	fwindow = quantum_hall(fn,r,n,x,y,f)
	fwindowp = [(fwindow[i], fwindow[i+1]) for i in range(0, len(fwindow), 2)]
	bnuwindow = [max(fwindowp[i])-min(fwindowp[i+1]) for i in range(0, len(fwindow)-1)]
	print(fwindow)
	print(fwindowp)
	plt.close('all')
	
	# Get Resistance
	e = 1.60217663*10**-19
	h = 6.62607015*10**-34
	Rk = h/(e**2)
	R = "Resistance"
	RName = "Klitzing "+ R 
	m = r[r["Measurement"]==n]
	I = m["I_ac"].values.tolist()[0]
	c = m["Correction"].values.tolist()[0]
	df = utilities.ReadDatFile(fn, correction=c, y=y)
	df[R] = (df[y].values/I)
	df[RName] = df[R].values/Rk
	
	# Make figure
	fig,ax = plt.subplots(2)
	
	ax[0].plot(df[x], df[R].values/1000, color='blue', label="Resistance")
	ax[0].set_ylabel(R + " (kOhms)")
	ax[1].plot(df[x], df[RName], color='blue', label="1/v")
	ax[1].set_ylabel(RName + " (Normalized)")
	
	#
	figxrange = max(df[x])-min(df[x])
	figxstep = figxrange/20
	Rfigrange = (max(df[R])-min(df[R]))/1000
	Vfigstep = Rfigrange/20
	
	Rfigrange2 = max(df[RName])-min(df[RName])
	Vfigstep2 = Rfigrange2/20
	
	# For every window that was elected
	for idx,i in enumerate(fwindowp):
		# Unpack window and get important values.
		even = idx%2
		if even:
			xloc = max(df[x])-5*figxstep
		else:
			xloc = min(df[x])+5*figxstep
		x1,x2 = i
		center = bnuwindow 	# Is this the Bnu value?
		cut = df[(df[x]>x1)&(df[x]<x2)]
		# Get the Nu value: mean nu over the cut window	
		v_value = np.mean(1/cut[RName].values)

		Ns = center*v_value*2.418*10**10
		
		# Get the Resistance in kOhms 
		yvals = cut[R].values/1000
		ymean = np.mean(yvals)	
		# Plot the kOhm markers.
		ax[0].axhline(y=ymean, color='red', linestyle='--')
		
		ax[0].text(xloc, ymean+Vfigstep/2, f"v: {v_value:0.1f} N_s: {Ns:0.2E}")
		
		# Get 1/v values. 
		yvals = cut[RName].values
		ymean = np.mean(yvals)
		# Plot the 1/v values.
		ax[1].axhline(y=ymean, color='red', linestyle='--')
		ax[1].text(xloc, ymean, f"v: {v_value:0.1f} N_s: {Ns:0.2E}")
	
	plt.show()

		

if __name__ == '__main__':
	f = [f'tristan_ringsgate-{i}.dat' for i in [17,14]]
