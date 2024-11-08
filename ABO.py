import pandas
import numpy as np
pandas.options.mode.chained_assignment = None 
import pylab as pl
from io import StringIO
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from scipy.optimize import curve_fit as fit
import matplotlib, peak_elector
import utilities, windower
#font = {'size'   : 8,
        #"family":'serif'}
#matplotlib.rc('font', **font)
#plt.rcParams.update({
        #"text.usetex": True})
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


fig_size_x, fig_size_y=8,8
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

def Extract(lst):
    return [item[0] for item in lst]

def ABOPeakElector(filename,		# String
                  df,			# Complete raw data DF
                  windows,		# list of tuples bounding coordinates, or 
                  record,		# Program-generated DF: Cols: window, rows: user x-click
                  bycoordinates,
                  x='B Field (T)',
                  y="P124A (V)",
                  selectregion=True
                  ):
    # Function Start
    try:
        for idx, w in enumerate(windows):
            print(idx,w[0])
            identifier=idx+w[0]
            if bycoordinates:
            	cut = df[(df[x]>w[0])&(df[x]<w[1])]
            else:
                cut = df.iloc[w[0]:w[1]]
            function="ThirdOrder"
            df, fig = GFF(cut, function, filename=f'{filename} cut {identifier}', window=w)
            if selectregion:
                print('\n\nClick the suspected peaks of the oscillation.')
                print('Press ENTER when all peaks have been denoted')
                print('Left click to mark. Right click to unmark.')
                tuples = plt.ginput(-1, timeout=False)
                plt.close()
            else:
                return fig,df
            td = Extract(tuples)
            record[idx] = td
            #record = pandas.concat([record,peak_elector.ElectPeaks(cut,\
			#f'{filename} cut {identifier}')], ignore_index=True, axis=1)
    except KeyboardInterrupt:
        with open("Record"+"".join(filename.split('.dat'))+".csv", 'w') as f:
            record.to_csv(f)
        raise KeyboardInterrupt("Interrupt recieved.")
    return record

def ABOMain(filenames:list, s,**kwargs):
    # Take a list of filenames, and a dict of correction factors, and do abo peak selection
    x=kwargs.get('x',"B Field (T)")
    results = kwargs.get('results',{})
    subwindowWidth=kwargs.get("subwindowWidth",None)
    numWindows=kwargs.get("numWindows",None)
    degenerate=kwargs.get("degenerate",False)
    # Make sure we were called correctly
    case1 = (numWindows is not None and subwindowWidth is not None)
    if subwindowWidth is None and numWindows is None:
        raise TypeError("Arguments: subwindowWidth and numWindows cannot BOTH be None")
    elif (case1 and degenerate==False) or (degenerate==True and subwindowWidth is None):
        raise TypeError("Arguments: numWindows cannot be passed with a\
			 subwindowWidth unless you want degenerate windowing. \
			 (Did you forget to pass degenerate=True?)")
    # Begin iteration
    for idx,filename in enumerate(filenames):
        # Search for old files.
        try:
             with open(f"ABO_Record{"".join(filename.split('.dat'))}.csv", 'r') as f:
                record = pandas.read_csv(f,index_col=0)
        except:
            record = pandas.DataFrame()
        # Open actual data file
        df = utilities.ReadDatFile(filename,correction=s[filename])
        if numWindows is None:
            # Make a coordinate-based window selection
            windows = windower.MakeSteppedSpan(df[x],stepsize=subwindowWidth)
            bycoordinates=True
        else:
            # Make a data-point based, index window selection
            splits, ranges = windower.SplitFile(df,numWindows,x=x)
            # Make another window?
            bycoordinates=False
        
        #
        try:
            record = ABOPeakElector(filename,df, windows,record,bycoordinates)
        except KeyboardInterrupt:
            with open("Record"+"".join(filename.split('.dat'))+".csv", 'w') as f:
                record.to_csv(f)
            raise KeyboardInterrupt()
        finally:
            with open("Record"+"".join(filename.split('.dat'))+".csv", 'w') as f:
                record.to_csv(f)

    return record

    
def main(): 
    k = [1,2,3,4,5,6,7,9]
    
    f = [f"tristan_ringsgate-{i}.dat" for i in k]
    correctionFactor = {}
    for i,v in enumerate(f):
        if k[i] in [3,4,5,6,7]:
            correctionFactor[v]=500/20*10**-6
        else:
            correctionFactor[v] = 1
    ABOMain(f,correctionFactor,subwindowWidth=4*5E-4)
    
if __name__ =='__main__':
    main()
