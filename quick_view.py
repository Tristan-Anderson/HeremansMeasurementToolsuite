import pandas as pd
import numpy as np
import pylab as pl
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
from utilities import FunctionInputHandler

def view(df, x='B Field (T)', y='P124A (V)'):
	""" 
	    Function which handles Hall Data. Requests singular file, with paramaters: x,y, and selx.
	    User can elect x-selection capabilities on call()
	"""
	# Force the user (future u) to call the function correctly
	df,f = FunctionInputHandler(df,f=df)
	
	# User selects data: show it to them.
	plt.plot(df[x], df[y], color='blue', label=y)
	
	plt.xlabel(x)
	plt.ylabel(y)
	
	# Titling.
	plt.title(f)
	plt.legend(loc='best', frameon=False)
	minx = min(df[x])
	miny=min(df[y])
	stepy = max(df[y])/miny
	plt.savefig(f'view_{f.split('.dat')[0]}')
	if __name__=='__main__':
		return False
	plt.close('all')

if __name__=='__main__':
	view("tristan_ringsgate-13.dat")
	plt.show()

