import pandas as pd
import numpy as np
import pylab as pl
import utilities, windower
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit

### General function to extract carrier densities from measurements
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
	print(m)
	I = m["I_ac"].values.tolist()[0]
	c = m["Correction"].values.tolist()[0]
	df = utilities.ReadDatFile(fn, correction=c, y=y)
	df[RName] = (df[y].values/I)/Rk
	dpts = []
	addition = [None,None]
	fig,ax = None,None
	while len(addition) ==2:
		addition = windower.selectWindower(df,x,\
				y,clickpoints=True, fig=fig,\
				ax=ax)
		dpts += addition
		fig,ax = plt.subplots(2)
		for x in dpts:
			ax[0].axvline(x=x, color='red', linestyle='--')
	cols = [i for i in range(len(dpts))]
	results = {i:dpts[i] for i in cols}
	rdf = pd.DataFrame(results)
	print(rdf)
	with open("QH_Record_"+fn.split('.dat')[0]+'.csv','w') as recordfile:
		rdf.to_csv(recordfile)		
	plt.close()

if __name__ == '__main__':
	f = [f'tristan_ringsgate-{i}.dat' for i in [8,11,14,17,18]]
	for i in f:
		print(i)
		with open(i, 'r') as fi:
			df =pd.read_csv(fi)
		hall(df,10**-6, selx=True)
		hall(df,10**-6,f=i, selx=True)
		hall(i, 10**-6, selx=True)
