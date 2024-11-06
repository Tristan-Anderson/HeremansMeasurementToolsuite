import numpy as np
import sys
import pandas as pd
from utilities import FunctionInputHandler


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
    assert max_val-min_val>=stepsize

    wholeIntervals= int(np.floor((max_val-min_val)/stepsize))-1
    if any([not isinstance(i, int) for i in iterable]):
    	spans = [[min_val+stepsize*i, min_val+stepsize*(i+1)-sys.float_info.epsilon] for i in range(wholeIntervals)]
    else:
    	spans = [[min_val+stepsize*i, min_val+stepsize*(i+1)-1] for i in range(wholeIntervals)]
    spans.append([min_val+stepsize*wholeIntervals,max_val])
    return spans

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

main()

def GetNumber(name):
    # Grabs a number from the filename that follows the lab's
    #   file naming convention
    a = name.split('-')
    b = a[-1]
    c = b.split('.')
    print(c[0], name)
    return c[0]

def ParseFile(file, s, y="P124A (V)"):
    # The correction factor is able to scale a specific
    #   stream of data by a constant. This is used ONLY
    #   when, during experiment, we don't enter the
    #   lock-in amplifier's sensitivity, properly into the
    #   LABView program.
    correction = s[int(GetNumber(file))]

    with open(file, 'r') as f:
        df = pandas.read_csv(f, delimiter=',')
    if isinstance(correction,float):
        df[y]=df[y].to_numpy()*correction

    return df
