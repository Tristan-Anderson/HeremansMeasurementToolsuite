import pandas as pd
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
from utilities import FunctionInputHandler
def line(m,x,b):
	return m*x+b

### General function to extract carrier densities from measurements
def hall(df, I, x='B Field (T)', y='P124A (V)', selx=False, f=None):
	""" 
	    Function which handles Hall Data. Requests singular file, with paramaters: x,y, and selx.
	    User can elect x-selection capabilities on call()
	"""
	# Force the user (future u) to call the function correctly
	df,f = FunctionInputHandler(df,f=f)
	
	# Default to fitting/plotting everything.
	#res = "Resistance (Ohm)"
	#df[res] = df[y].values / I
	datax = df[x]
	datay = df[y]
	xd = df[x].values
	yd = df[y].values
	
	if selx:
		# User selects data: show it to them.
		plt.scatter(df[x], df[y], color='blue')
		print("Delinate fit region between two clicks.")
		uinput = plt.ginput(2, timeout=0)
		# Only 2 allowed.
		plt.close('all')
		if len(uinput) == 2:
			# Extract x's
			xs = sorted([float(i[0]) for i in uinput])

			# Create !df[x] for ploting 
			dfplt = df[~(df[x]>xs[0])|~(df[x]<xs[1])]
			datax, datay = dfplt[x], dfplt[y]
			
			# Create df for fit.
			dff = df[(df[x]>xs[0])&(df[x]<xs[1])]
			xd, yd = dff[x], dff[y]
			
			# Flash user
			plt.scatter(xd,yd, color='blue')
			plt.scatter(datax,datay, color='black')
			plt.pause(.75)
			plt.close('all')
		else:
			pass
	
	# Fit data
	pvar,_ = fit(line,xd,yd)
	
	# Grab original data to plot fit-line with.
	xfit = df[x].values
	xfit = np.linspace(min(df[x].values), max(df[x].values), num=500)
	yfit = line(xfit, *pvar)
	
	# Plot data you've gathered	
	plt.scatter(datax,datay, color='black', label='dataset', s=3)
	plt.scatter(xd,yd, color='blue', label='data for hall fit', s=3)
	plt.plot(xfit, yfit, color='red', label='fit')
	ytot = list(yd)+list(yfit)
	xtot = list(xfit)+list(xd)
	
	# Calculate density assuming e
	Ns = I/(-1.602*10**(-19) *pvar[0])
	print(f'{f} e- density: {Ns:0.4E}')
	plt.text(min(xtot), min(ytot), r"N_s: " f"{Ns:0.3E}"+ r"$\frac{1}{m^2}$")
	plt.xlabel(x)
	plt.ylabel(y)
	
	# Titling.
	plt.title(f)
	plt.legend(loc='best', frameon=False)
	print(f)
	plt.savefig(f'hall_{f.split('.dat')[0]}')
	plt.close('all')


if __name__ == '__main__':
	f = [f'tristan_ringsgate-{i}.dat' for i in [8,11,14,17,18]]
	for i in f:
		print(i)
		with open(i, 'r') as fi:
			df =pd.read_csv(fi)
		hall(df,10**-6, selx=True)
		hall(df,10**-6,f=i, selx=True)
		hall(i, 10**-6, selx=True)
