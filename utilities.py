import pandas as pd

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
