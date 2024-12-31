import pandas
import pandas as pd
import numpy as np
import windower as _windower
from matplotlib import pyplot as plt

fig_size_x, fig_size_y = 5, 8
def nearest(iterable, test):
	pass

def Linear(x,m,b):
    return m*x+b
def SecondOrder(x,a,b,c):
    return x*Linear(x,a,b)+c
def ThirdOrder(x,a,b,c,d):
    return x*SecondOrder(x,a,b,c)+d
def Sin(x,a,p,w,y0):
    return a*np.sin(w*x-p)+y0

def GFF(df, function, **kwargs):
    #  https://github.com/Tristan-Anderson/Slifer_Lab_NMR_Toolsuite in NMR_Analyzer.py
    """
    Generalized Fitting Function


    """
    def get_function(f_name,xdata,var):
        if f_name == "Sin":
            yfit = Sin(xdata, var[0], var[1], var[2], var[3])

        elif f_name == "ThirdOrder":
            yfit = ThirdOrder(xdata, var[0], var[1], var[2], var[3])

        elif f_name == "SecondOrder":
            yfit = SecondOrder(
                xdata, var[0], var[1], var[2])
        elif f_name=="None":
            yfit = df[y].values

        return yfit

    # Get kw arguments
    window = kwargs.pop("window",[0,0])
    # Sets which Y to evaluate
    y = kwargs.pop('y', "P124A (V)")

    # Sets which X to evaluage
    x = kwargs.pop('x', "B Field (T)")

    # Toggles saving the figure
    savefit = kwargs.pop('savefit', False)

    # If the above is toggled, the file needs a name
    filename = kwargs.pop('filename', 'UNNAMED_GRAPH')

    #automated?
    automated = kwargs.pop('automated', False)

    # selectregion?
    selectregion = kwargs.pop('selectregion', False)


    # Sets the window for data
    xmin = kwargs.pop('xmin', None)
    xmax = kwargs.pop('xmax', None)
    ogdf = df.copy(deep=True)
    # This requires an index
    if any(isinstance(i, float) for i in window):
        df = df[(df[x]>window[0]) & (df[x]<window[1])]
    else:
    	df=df.iloc[window[0]:window[1]]
    cut_data = pandas.DataFrame()
    if selectregion:
        fig, ax = plt.subplots(figsize=(fig_size_x, fig_size_y))
        ax.scatter(df[x],df[y],label=y,color='blue')
        ax.legend(loc='best')
        tuples = plt.ginput(2,timeout=False)
        plt.close()
        xlims = numpy.asarray([c[0] for c in tuples])
        cut_data = df[(df[x] > xlims[0]) & (df[x] < xlims[1])]

    elif xmin is not None and xmax is not None:
        cut_data = df[(df.x > xmin) & (df.x < xmax)]

    else:
        fit_data = df

    #print(df)

    fit_data = df.drop(cut_data.index.to_numpy())

    #print(fit_data)


    xdata = fit_data[x].values
    ydata = fit_data[y].values


    var, _ = fit(eval(function), xdata, ydata)
    # Fit that stuff
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"

    df.insert(column=fitname, value=get_function(function, df[x], var),loc=4)
    df.insert(column=fitsub, value=df[y]-df[fitname],loc=5)

    if not automated:

        fig, ax = plt.subplots(3,figsize=(fig_size_x, fig_size_y))
        #plt.title(filename)
        ax[0].set_title(filename)
        ax[0].scatter(ogdf[x],ogdf[y],label='data',color='blue')
        ax[0].scatter(df[x],df[y],label='current selection',color='red')

        ax[1].plot(df[x],df[fitname],label='fit',color="red")
        ax[1].plot(df[x],df[y],label='data',color='blue')

        ax[2].plot(df[x],df[fitsub],label='fit subtraction',color="black")

        for i in ax:
            i.grid(True)
            i.legend(loc="best")

    return df, fig

def generateWindowRecord(windows,filename,prefix):
	d = {}
	d['window num']=[]
	d['x1']=[]
	d['x2']=[]
	for i,(x1,x2) in enumerate(windows):
		d['x1'].append(x1)
		d['x2'].append(x2)
		d['window num'].append(i)
	df = pd.DataFrame(d)
	with open(f"WINDOWS_{prefix}_{"".join(filename.split('.dat'))}.csv",'w') as f:
		df.to_csv(f,index=False)

def readWindows(filename,prefix):
	with open(f"WINDOWS_{prefix}_{"".join(filename.split('.dat'))}.csv",'r') as f:
		df= pd.read_csv(f)
	return df
		

def getwindows(filename,s,prefix,**kwargs):
    x=kwargs.get('x',"B Field (T)")
    results = kwargs.get('results',{})
    subwindowWidth=kwargs.get("subwindowWidth",None)
    numWindows=kwargs.get("numWindows",None)
    degenerate=kwargs.get("degenerate",False)
    # Search for old files.
    try:
         with open(f"{prefix}_Record{"".join(filename.split('.dat'))}.csv", 'r') as f:
            record = pandas.read_csv(f,index_col=0)
    except:
        record = pandas.DataFrame()
    # Open actual data file
    df = ReadDatFile(filename,correction=s[filename])
    try:
        if numWindows is None:
            # Make a coordinate-based window selection
            windows = _windower.MakeSteppedSpan(df[x],stepsize=subwindowWidth)
            bycoordinates=True
        else:
            # Make a data-point based, index window selection
            splits, ranges = _windower.SplitFile(df,numWindows,x=x)
            # Make another window?
            bycoordinates=False
    except AssertionError:
        print(f"Assertion error detected. Meaning df.x.max - df.x.min < stepsize! Skipping!")
        print(f"{filename} has xmin: {min(df[x]):0.1E}, xmax: {max(df[x]):0.1E}, range: {max(df[x])-min(df[x]):0.1E}, stepsize: {subwindowWidth:0.1E} (T)")
    generateWindowRecord(windows,filename,prefix)
    return windows, bycoordinates, df, record
	
def FunctionInputHandler(df,f=None):
	# Force the user (future u) to call functions correctly
	if isinstance(df, pd.DataFrame):
		if f is None:
			f = str(input("Enter Measurement name: "))
		elif isinstance(f, str):
			pass
		else:
			raise TypeError("Dataframe was passed, but no original file name.")
	elif isinstance(df, str):
		f = df
		with open(f, 'r') as _f:
			df = pd.read_csv(_f)
	else:
		raise TypeError("called function expects a string path to a csv, or dataframe.")
	return df, f


def GetNumber(name):
    # Grabs a number from the filename that follows the lab's
    #   file naming convention: {NAME}-{NUM}.dat
    a = name.split('-')
    b = a[-1]
    c = b.split('.')
    print(c[0], name)
    return c[0]


def ReadDatFile(file, correction=1, y="P124A (V)"):
    """
        Reads the file, corrects data if necessary
	c: correction factor to apply on y-column
    """
    # The correction factor is able to scale a specific
    #   stream of data by a constant. This is used ONLY
    #   when, during experiment, we don't enter the
    #   lock-in amplifier's sensitivity, properly into the
    #   LABView program.
    df,f = FunctionInputHandler(file)

    df[y]=df[y].to_numpy()*correction

    return df

def Sin(x, a,b,c,d):
    return a*numpy.sin(b*x-c)+d

def ThirdOrder(x, a,b,c,d):
    return a*x**3+b*x**2+c*x+d

def SecondOrder(x, a,b,c):
    return a*x**2+b*x+c


def GFF(df, function, **kwargs):
    #  https://github.com/Tristan-Anderson/Slifer_Lab_NMR_Toolsuite in NMR_Analyzer.py
    """
    Generalized Fitting Function


    """
    def get_function(f_name,xdata,var):
        if f_name == "Sin":
            yfit = Sin(xdata, var[0], var[1], var[2], var[3])

        elif f_name == "ThirdOrder":
            yfit = ThirdOrder(xdata, var[0], var[1], var[2], var[3])

        elif f_name == "SecondOrder":
            yfit = SecondOrder(
                xdata, var[0], var[1], var[2])

        return yfit

    # Get kw arguments
    window = kwargs.pop("window",[0,0])
    # Sets which Y to evaluate
    y = kwargs.pop('y', "P124A (V)")

    # Sets which X to evaluage
    x = kwargs.pop('x', "B Field (T)")

    # Toggles saving the figure
    savefit = kwargs.pop('savefit', False)

    # If the above is toggled, the file needs a name
    filename = kwargs.pop('filename', 'UNNAMED_GRAPH')

    #automated?
    automated = kwargs.pop('automated', False)

    # selectregion?
    selectregion = kwargs.pop('selectregion', False)


    # Sets the window for data
    xmin = kwargs.pop('xmin', None)
    xmax = kwargs.pop('xmax', None)
    ogdf = df.copy(deep=True)
    df=df.iloc[window[0]:window[1]]
    cut_data = pandas.DataFrame()
    if selectregion:
        fig, ax = plt.subplots(figsize=(fig_size_x, fig_size_y))
        ax.scatter(df[x],df[y],label=y,color='blue')
        ax.legend(loc='best')
        tuples = plt.ginput(2,timeout=False)
        plt.close()
        xlims = numpy.asarray([c[0] for c in tuples])
        cut_data = df[(df[x] > xlims[0]) & (df[x] < xlims[1])]

    elif xmin is not None and xmax is not None:
        cut_data = df[(df.x > xmin) & (df.x < xmax)]

    else:
        fit_data = df

    #print(df)

    fit_data = df.drop(cut_data.index.to_numpy())

    #print(fit_data)


    xdata = fit_data[x].values
    ydata = fit_data[y].values


    var, _ = fit(eval(function), xdata, ydata)
    # Fit that stuff
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"

    df.insert(column=fitname, value=get_function(function, df[x], var),loc=4)
    df.insert(column=fitsub, value=df[y]-df[fitname],loc=5)

    if not automated:

        fig, ax = plt.subplots(3,figsize=(fig_size_x, fig_size_y))
        plt.title(filename)

        ax[0].scatter(ogdf[x],ogdf[y],label='data',color='blue')
        ax[0].scatter(df[x],df[y],label='current selection',color='red')

        ax[1].plot(df[x],df[fitname],label='fit',color="red")
        ax[1].scatter(df[x],df[y],label='data',color='blue')

        ax[2].scatter(df[x],df[fitsub],label='fit subtraction',color="black")

        for i in ax:
            i.grid(True)
            i.legend(loc="best")

    return df, fig

def pairwise(a):
    # len(a) choose 2
    b = len(a)
    res = []
    for i in range(b):
        sel = a[i]
        left = a[i+1:]
        if len(left) == 0:
            break
        tl = [(sel,l) for l in left]
        res += tl
    return res

def quadwise(contacts):
	from pprint import pprint
	# len(contacts) choose 4
	# len_c choose 4 = (len_c choose 2)(len_c-2 choose 2)

	# Establish V pairs (aka all 2 ports)
	pair1 = pairwise(contacts)
	final = []
	
	# For each tuple formed by (N choose 2)
	for (x1,x2) in pair1:
		tl = []
		candidates = []
		# Find the contacts that are not in in the (N choose 2)
		for c in contacts:
			if c not in [x1,x2]:
				candidates.append(c)
		# Then formulate (N choose 2)'s compliment
		# Being (N-2 choose 2)
		n_2c2 = pairwise(candidates)
		tl = [(x1,x2,x3,x4) for (x3,x4) in n_2c2]
		final += tl
	print("2-port combinations", len(pair1))
	pprint(pair1)
	print("4-port combinations",len(final))
	pprint(final)
