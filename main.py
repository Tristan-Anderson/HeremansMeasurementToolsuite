import pandas, numpy
pandas.options.mode.chained_assignment = None 
import pylab as pl
from io import StringIO
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from scipy.optimize import curve_fit as fit
import matplotlib
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

def Sin(x, a,b,c,d):
    return a*numpy.sin(b*x-c)+d

def ThirdOrder(x, a,b,c,d):
    return a*x**3+b*x**2+c*x+d

def SecondOrder(x, a,b,c):
    return a*x**2+b*x+c


def GetNumber(name):
    # Grabs a number from the filename that follows the lab's
    #   file naming convention
    a = name.split('-')
    b = a[-1]
    c = b.split('.')
    print(c[0], name)
    return c[0]

def ParseFile(file, s, y="P124A (V)"):
    """
        Reads the file, corrects data if necessary (sensitivity scales everything such 
        that the final reported value is in Volts.)
    """
    # The correction factor is able to scale a specific
    #   stream of data by a constant. This is used ONLY
    #   when, during experiment, we don't enter the 
    #   lock-in amplifier's sensitivity, properly into the
    #   LABView program.
    correction = s[int(GetNumber(file))]
    
    with open(file, 'r') as f:
        df = pandas.read_csv(f, delimiter=',')
    if type(correction) ==float:
        df[y]=df[y].to_numpy()*correction

    return df

def SplitFile(*args,**kwargs):
    df= args[0]
    numWindows=args[1]

    B=kwargs.get('B',"B Field (T)")
    width=kwargs.get("width", int(len(df[B])/numWindows))
    degenerate=kwargs.get('degenerate',False)
    if degenerate:
        stepsize=kwargs.get('stepsize',1)

    # Define split-points for the DF based on the index
    #   Split points correspond to the nth datapoint.
    splits_for_df = numpy.arange(0,len(df[B])-1, width)
    
    # range_for_df is a nested-list [[start1, end1], ...]
    #   the difference between endX and start(X+1) is 1.
    #   No data is "sliced out" of analysis.
    
    if degenerate:
        range_for_df = MakeSteppedSpan(splits_for_df, len(df[B])-1, 50,degenerate=degenerate,width=width)
    else:
        range_for_df = MakeSplitSpan(splits_for_df,max_val=len(df[B])-1)

    return splits_for_df,range_for_df

def USelect(df,T,B,y,file,s):
    b=""
    while b !='y':
        # Untill the user is satisfied, allow them to select
        #       Data that they need.
        plt.close('all')
        plt.plot(df[B], df[y], label="Data")
        plt.grid(True)
        plt.legend(loc='best')
    
        coords = plt.ginput(2)
        xax = [c[0] for c in coords]
        cut = df[(df[B]<max(xax)) & (df[B]>min(xax))]
        other = df[df[B]>max(xax)]
        other = pandas.concat([other, df[df[B]<min(xax)]])
        print(other)
        plt.close('all')
        fig,ax = plt.subplots(2)
        ax[0].scatter(other[B], other[y])
        ax[1].scatter(cut[B], cut[y])
        for i in ax:
            i.grid(True)
            i.legend(loc='best')
        plt.show()
        plt.close('all')
        print("Window Width: ",max(xax)-min(xax))
        b = input("Ok? (y/n): ")
    Analyze(df,0,T,B,y,file,cut,other,s)

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

def Extract(lst):
    return [item[0] for item in lst]

    
def ABOscillation(filenames,s,windows,df,record,succeded, **kwargs):
    try:
        for idx, w in enumerate(windows):
            identifier=idx+w[0]
            record = pandas.concat([record,ElectPeaks(df, w,  i,identifier)], ignore_index=True, axis=1)
            print(record)
    except KeyboardInterrupt:
        with open("Record"+"".join(filenames[0].split('.dat'))+".csv", 'w') as f:
            record.to_csv(f)
    return record

def thirdOrderFitSub(df,**kwargs):
    y = kwargs.pop('y', "P124A (V)")
    x = kwargs.pop('x', "B Field (T)")
    function="ThirdOrder"
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"
    X,Y = df[x],df[y]
    pvar,_=fit(ThirdOrder,X,Y)
    df[fitname] = ThirdOrder(X,*pvar)
    df[fitsub] = df[y]-df[fitname]
    return df

def TileResultsGrabber(filename,s,windows,df,record,**kwargs):
    y = kwargs.pop('y', "P124A (V)")
    x = kwargs.pop('x', "B Field (T)")
    results = kwargs.pop('results', {})
    function="ThirdOrder"
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"
    results = kwargs.get('results',{})
    for idx,w in enumerate(windows):
        selected=record[str(idx)].values # User elected peak locations
        #selected=record[idx].values # User elected peak locations
        selected = selected[~numpy.isnan(selected)]
        print(filename)
        if selected.size == 0:
            print("Window", idx, "Empty. Skipping.")
            continue
        cut = df.iloc[w[0]:w[1]]
        cut = thirdOrderFitSub(cut)
        ydata = cut[fitsub].values
        xdata = cut[x].values
        results["".join(filename.split('.dat'))+'-'+str(idx)] = {fitsub:ydata, x:xdata,'election':selected}
    return results

def Toggle(filenames,s,**kwargs):
    B=kwargs.get('B',"B Field (T)")
    replot=kwargs.get('replot',False)
    oscillation=kwargs.get('oscillation',False)
    subwindowWidth=kwargs.get("subwindowWidth",2/5*1E-2)
    numWindows=kwargs.get("numWindows",10)
    results = kwargs.get('results',{})
    succeded = 0
    print(filenames, filenames[-1], 366)
    for idx,i in enumerate(filenames):
        try:
             with open("Record"+"".join(i.split('.dat'))+".csv", 'r') as f:
                record = pandas.read_csv(f,index_col=0)
        except:
            record = pandas.DataFrame()
        df = ParseFile(i,s)
        # No overlap in splits.
        BM,Bm= max(df[B]),min(df[B])
        numWindows=abs(numpy.ceil((BM-Bm)/subwindowWidth))
        # window = [[start index df window:int, end index df window:int], ... ]
        splitpoints, windows = SplitFile(df,numWindows)
        print("window indecies:",windows)
        xax=df[B]
        print("points per subwindow",int(len(xax)/numWindows))
        # Identifier need to uniquely track the nth subwindow across all windows.
        if replot:
            succeded += ReplotData(filenames,s,windows, df,record,succeded,**kwargs)
        elif oscillation:
            record = ABOscillation(filenames,s,windows,df,record,succeded,**kwargs)
            with open("Record"+"".join(i.split('.dat'))+".csv", 'w') as f:
                record.to_csv(f)
        else:
            results = TileResultsGrabber(i,s,windows,df,record,results=results)
    if replot:
        return succeded
    elif oscillation:
        return record
    else:
        return results

def ReplotData(filenames,s,windows, df,record,succeded,**kwargs):
    B=kwargs.get('B',"B Field (T)")
    subwindowWidth=kwargs.get("subwindowWidth",10E-4)
    numWindows=kwargs.get("numWindows",10)
    for idx, w in enumerate(windows):
        cut_df = df.iloc[w[0]:w[1]]
        identifier=idx+w[0]
        fig,_ = ElectPeaks(df, w,  i, identifier,selectregion=False)
        fig.tight_layout(pad=1.75)
        ax = fig.axes
        ax[0].title.set_text(i)
        ax[1].title.set_text("Window "+str(idx)+" Raw Data and Fit")
        ax[2].title.set_text("Fit subtracted window with user-elected peaks")
        selected=record[str(idx)].values # User elected peak locations
        if selected.size == 0:
            print("Window", idx, "Empty. Skipping.")
        else:
            for xt in selected:
                ax[2].axvline(x=xt,color='red',linestyle='--') # Marking peaks.
            for z in ax:
                z.set_ylabel("Volts (V)")
                z.set_xlabel("Magnetic Field (T)")
                z.grid(False)
            #plt.show()
            #fig.set_size_inches(6.5,9)
            #plt.savefig("".join(filenames[0].split('.dat'))+"_Window_"+str(idx),dpi=150)
            #plt.show()
            plt.close('all')
            succeded+=1
    return succeded

def TileResults(r,**kwargs):
    y = kwargs.pop('y', "P124A (V)")
    x = kwargs.pop('x', "B Field (T)")
    function="ThirdOrder"
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"
    print(len(r))
    plt.rcParams["figure.figsize"] = (18,32)
    ncols,nrows =3,12
    fig,ax = plt.subplots(ncols=ncols,nrows=nrows)
    l=0
    m=0
    f =0
    fig2, ax2 = plt.subplots(1,figsize=(10,5))
    for i,k in enumerate(r.keys()):
        if i % ncols==0 and i!=0:
            l=0
            m+=1
        print(m,l)
        print(k)
        v = r[k]
        Xdata=v[x]
        Ydata=v[fitsub]/(100*10**-9)
        signs = Ydata/abs(Ydata)
        election=v['election']
        if f==12:
            ax2.scatter([k*1000 for k in Xdata],Ydata,color='black')
            for xt in election:
                ax2.axvline(x=xt*1000,color='red',linestyle='--') # Marking peaks.
            ax2.set_xlabel("B (mT)")
            ax2.set_ylabel("R (Ω)")
            plt.show()
            fig2.savefig("SelectedOscillation.svg")
            fig2.savefig("SelectedOscillation")
            exit()
        ax[m,l].scatter(Xdata,Ydata,color='black')
        for xt in election:
                ax[m,l].axvline(x=xt,color='red',linestyle='--') # Marking peaks.
        ax[m,l].grid(True)
        l+=1
        f+=1
    
            
    fig.suptitle("Regions of Measurements Where Oscillations were Found")
    #plt.show()
    plt.savefig("Tiled.svg",dpi=180)
    
    
    

k = [1,2,3,4,5,6,7,9]
sensitivity = {i:[1] for i in range(12)}
for i in range(1,3):
	sensitivity[i]=500*10**-6
for i in range(3,7):
	sensitivity[i]=500/20*10**-6
for i in range(7,9):
	sensitivity[i]=100*10**-6
sensitivity[9]=200*10**-6

f = [f"tristan_ringsgate-{i}.dat" for i in k]
j = 0
results = {}
Toggle([f[-1]],sensitivity,oscillation=True)
"""
for i in f:
    j+=Toggle([i],sensitivity,oscillation=True)
    results = results | Toggle([i],sensitivity,results=results)
#print(results)
#TileResults(results)
"""
