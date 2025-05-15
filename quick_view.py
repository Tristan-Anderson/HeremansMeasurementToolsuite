import pylab as pl
from matplotlib import pyplot as plt
from utilities import FunctionInputHandler

def view(df, x='B Field (T)', y='P124A (V)'):
	""" 
	    Function which handles Hall Data. Requests singular file, with paramaters: x,y, and selx.
	    User can elect x-selection capabilities on call()
	"""
	# Force the user (future u) to call the function correctly
	df,f = FunctionInputHandler(df,f=df)
	pl.rcParams['figure.figsize']=8,8
	
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
	plt.savefig(f'view_{f.split('.dat')[0]}',dpi=300)
	if __name__=='__main__':
		return False
	plt.close('all')

def electView(df,xi,xf, x='B Field (T)', y='P124A (V)'):
	df,f = FunctionInputHandler(df,f=df)
	df = df[(df[x]>xi) & (df[x]<xf)]

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
	view("tristan_ringsgate-13.dat")
	plt.show()

