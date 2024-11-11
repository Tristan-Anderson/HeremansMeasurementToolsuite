import numpy as np
import sys
import pandas as pd
from matplotlib import pyplot as plt
from utilities import FunctionInputHandler

def selectWindower(df,x,y, clickpoints=False, fig=None, ax=None):
    b=""
    cpts = []
    ifg, iax = fig, ax
    # Untill the user is satisfied, allow them to select
    #       Data that they need.
    if fig is None:
        plt.close('all')
        fig,ax = plt.subplots(1)
        ax.plot(df[x], df[y], label="Data")
        ax.grid(True)
        plt.legend(loc='best', frameon=False)
    else:
        plt.close('all')
        ax[0].scatter(df[x], df[y],color='blue')
        ax[1].scatter(df[x], df[y],color='blue')
    coords = plt.ginput(2,timeout=0)
    if len(coords) == 0:
        return []
    xax = [c[0] for c in coords]
    cut = df[(df[x]<max(xax)) & (df[x]>min(xax))]
    other = df[~(df[x]<max(xax)) | ~(df[x]>min(xax))]
    plt.close('all')
    fig,ax = plt.subplots(2)
    ax[0].scatter(other[x], other[y],color='blue')
    ax[0].scatter(cut[x], cut[y],color='red')
    ax[1].scatter(cut[x], cut[y],color='red')
    for i in ax:
        i.grid(True)
        i.set_xlabel(x)
        i.set_ylabel(y+" (1/v)")
    while True:
        c2 = plt.ginput(2,timeout=0)
        if len(c2) == 0:
            return []
        else:
            plt.close('all')
            fig,ax = plt.subplots(2)
            xax = [c[0] for c in c2]
            cut = df[(df[x]<max(xax)) & (df[x]>min(xax))]
            other = df[~(df[x]<max(xax)) | ~(df[x]>min(xax))]
            ax[0].scatter(other[x], other[y],color='blue')
            ax[0].scatter(cut[x], cut[y],color='red')
            ax[1].scatter(cut[x], cut[y],color='red')
        for i in ax:
            i.grid(True)
            i.set_xlabel(x)
            i.set_ylabel(y+" (1/v)")
            
    if clickpoints:
        return xax
    else:
        return cut

def MakeSteppedSpan(iterable, **kwargs):
    """
        take a minimum and maximum value from a list,
        Use their difference for a width.
        From minimum to maximum value, create tuples that
        span width from minimum to maximum value.
	For ints: { [i_0, i_1], [i_1+1, i_2] ... }
	For floats: { [i_0, i_1-epsilon], [i_1, i_2] }
        return tuples
    """
    min_val = kwargs.pop("min_val", min(iterable))
    max_val = kwargs.pop("max_val", max(iterable))
    stepsize = kwargs.pop("stepsize",0)
    assert abs(max_val-min_val)>=stepsize
    epsilon = sys.float_info.epsilon

    wholeIntervals= int(np.floor((max_val-min_val)/stepsize))-1
    if any([not isinstance(i, int) for i in iterable]):
    	spans = [[min_val+stepsize*i, min_val+stepsize*(i+1)\
                  -epsilon] for i in range(wholeIntervals)]
    else:
    	spans = [[min_val+stepsize*i,\
                   min_val+stepsize*(i+1)-1] for i in range(wholeIntervals)]
    spans.append([min_val+stepsize*wholeIntervals,max_val])
    return spans

def SplitFile(df, numWindows, x="B Field (T)",**kwargs):
    width=kwargs.get("width", int(len(df[B])/numWindows))
    degenerate=kwargs.get('degenerate',False)
    if degenerate:
        try:
        	stepsize=kwargs.get('stepsize')
        except Exception as e:
                raise KeyError("degenerate kwarg must be passed with stepsize:int kwarg")

    # Define split-points for the DF based on the index
    #   Split points correspond to the nth datapoint.
    splits_for_df = numpy.arange(0,len(df[x])-1, width)

    # range_for_df is a nested-list [[start1, end1], ...]
    #   the difference between endX and start(X+1) is 1.
    #   No data is "sliced out" of analysis.

    if degenerate:
        range_for_df = MakeSteppedSpan(splits_for_df, len(df[x])-1,\
			      50,degenerate=degenerate,width=width)
    else:
        range_for_df = MakeSplitSpan(splits_for_df,max_val=len(df[x])-1)

    return splits_for_df,range_for_df

def main():
	x='B Field (T)'
	with open("tristan_ringsgate-1.dat", 'r') as f:
		df = pd.read_csv(f)
	x = df[x].values
	ss = 1E-2
	print(max(x), ss)
	span = MakeSteppedSpan(x,stepsize=ss)
	print('nums:', span)
	xidx = [i for i in range(len(x))]
	print(max(xidx), 25)
	span = MakeSteppedSpan(xidx,stepsize=25)
	print('ints:', span)

def GetNumber(name):
    # Grabs a number from the filename that follows the lab's
    #   file naming convention
    a = name.split('-')
    b = a[-1]
    c = b.split('.')
    print(c[0], name)
    return c[0]


if __name__ == '__main__':
	main()
