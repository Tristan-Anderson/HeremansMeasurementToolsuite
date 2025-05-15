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
    #xfit = np.linspace(min(df[x].values), max(df[x].values), num=500)
    yfit = line(xfit, *pvar)

    # Plot data you've gathered
    plt.scatter(datax,datay*10**3, color='black', label='dataset', s=3)
    plt.scatter(xd,yd*10**3, color='blue', label='data for hall fit', s=3)
    plt.plot(xfit, yfit*10**3, color='red', label='fit')
    ytot = list(yd)+list(yfit)
    xtot = list(xfit)+list(xd)

    # Calculate density assuming e
    Ns = I/(-1.602*10**(-19) *pvar[0])
    print(f'{f} e- density: {Ns:0.4E}')
    print(min(xtot), min(ytot))
    plt.text(min(xtot), min(ytot)*10**3, r"$N_s$: " f"{Ns:0.3E}"+ r"$\frac{1}{m^2}$")
    plt.xlabel(x)
    plt.ylabel(y + "(mV)")

    # Titling.
    plt.title(f)
    plt.legend(loc='best', frameon=False)
    print(f)
    plt.savefig(f'hall_{f.split('.dat')[0]}')
    plt.close('all')


if __name__ == '__main__':
    #f = [f'tristan_ringsgate-{i}.dat' for i in [8,11,14,17,18]]
    import pylab as pl

    pl.rcParams['figure.figsize'] = 8,8
    pl.rcParams['lines.linewidth'] = 1.5
    pl.rcParams["figure.autolayout"] = True
    pl.rcParams['font.family'] = 'serif'
    pl.rcParams['font.weight'] = 'normal'
    pl.rcParams['font.size'] = 12

    pl.rcParams['font.sans-serif'] = 'serif'
    pl.rcParams['text.usetex'] = False
    pl.rcParams['axes.linewidth'] = 1.5
    pl.rcParams['axes.titlesize'] = 'medium'
    pl.rcParams['axes.labelsize'] = 'medium'

    pl.rcParams['xtick.major.size'] = 8
    pl.rcParams['xtick.minor.size'] = 4
    pl.rcParams['xtick.major.pad'] = 8
    pl.rcParams['xtick.minor.pad'] = 8
    pl.rcParams['xtick.color'] = 'k'
    pl.rcParams['xtick.labelsize'] = 'small'
    pl.rcParams['xtick.direction'] = 'in'

    pl.rcParams['ytick.major.size'] = 8
    pl.rcParams['ytick.minor.size'] = 4
    pl.rcParams['ytick.major.pad'] = 8
    pl.rcParams['ytick.minor.pad'] = 8
    pl.rcParams['ytick.color'] = 'k'
    pl.rcParams['ytick.labelsize'] = 'medium'
    pl.rcParams['ytick.direction'] = 'in'
    path = r"C:\Users\Tristan\Desktop\work\Research\heremans\MonarkMeasurements\20250203Cooldown1AmpOhmn\Hall"
    f = [f'{path}\\combined.csv']
    for i in f:
        print(i)
        with open(i, 'r') as fi:
            df =pd.read_csv(fi)
        df = df.dropna()
        df["B (T)"] = df["AMI_430_Z_Axis_QCoDeS - B"]
        df["Resistance"] = df["lockin_1 - X"]/(100*10**-9)
        hall(df,100*10**-9, selx=True,x="AMI_430_Z_Axis_QCoDeS - B", y="lockin_1 - X")