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
        selected = selected[~np.isnan(selected)]
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

def ABOPeakElector(filename,df,windows,record,bycoordinates):
    try:
        for idx, w in enumerate(windows):
            identifier=idx+w[0]
            if bycoordinates:
            	cut = df[(df[x]>w[0])&(df[x]<w[1])]
            else:
                cut = df.iloc[w[0]:w[1]]
            record = pandas.concat([record,peak_elector.ElectPeaks(cut)], ignore_index=True, axis=1)
    except KeyboardInterrupt:
        with open("Record"+"".join(filename.split('.dat'))+".csv", 'w') as f:
            record.to_csv(f)
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
        raise TypeError("Arguments: numWindows cannot be passed with a subwindowWidth unless you want degenerate windowing. (Did you forget to pass degenerate=True?)")
    # Begin iteration
    for idx,filename in enumerate(filenames):
        # Search for old files.
        try:
             with open(f"ABO_Record{"".join(filename.split('.dat'))}.csv", 'r') as f:
                record = pandas.read_csv(f,index_col=0)
        except:
            record = pandas.DataFrame()
        # Open actual data file
        print(s)
        df = utilities.ReadDatFile(filename,correction=s[filename])
        if numWindows is None:
            # Make a coordinate-based window selection
            windows = windower.MakeSteppedSpan(df[x].values,stepsize=subwindowWidth)
            bycoordinates=True
        else:
            # Make a data-point based, index window selection
            splits, ranges = windower.SplitFile(df,numWindows,x=x)
            bycoordinates=False
        
        #
        record = ABOPeakElector(filename, windows,df,record,bycoordinates)
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
    ABOMain(f,correctionFactor,subwindowWidth=4*5E-5)
    
if __name__ =='__main__':
    main()
