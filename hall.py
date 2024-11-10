import pandas as pd
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
from utilities import FunctionInputHandler
pl.rcParams['figure.figsize']  = 8.5, 11
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'normal'
pl.rcParams['font.size']       = 12

pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = False
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad']  = 8
pl.rcParams['xtick.minor.pad']  = 8
pl.rcParams['xtick.color']      = 'k'
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad']  = 8
pl.rcParams['ytick.minor.pad']  = 8
pl.rcParams['ytick.color']      = 'k'
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'
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
	
	# Calculate density assuming e
	Ns = I/(-1.602*10**(-19) *pvar[0])
	print(f'{f} e- density: {Ns:0.4E}')
	
	# Titling.
	plt.title(f)
	plt.legend(loc='best', frameon=False)
	print(f)
	plt.savefig(f'hall_{f.split('.dat')[0]}')
	plt.close('all')


print(__name__)
if __name__ == '__main__':
	f = [f'tristan_ringsgate-{i}.dat' for i in [8,11,14,17,18]]
	for i in f:
		print(i)
		with open(i, 'r') as fi:
			df =pd.read_csv(fi)
		hall(df,10**-6, selx=True)
		hall(df,10**-6,f=i, selx=True)
		hall(i, 10**-6, selx=True)
