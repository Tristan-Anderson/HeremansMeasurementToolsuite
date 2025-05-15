import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#Temp (K),B Field (T),P124A (V),SR830 X (V),SR830 Y (V),EGG7265 X (V),EGG7265 Y (V),K6487 V (V),K6487 I (A)


def ElectPeaks(df, filename, **kwargs):
    # Get column labels (x,y), assign values to X, Y
    y = kwargs.pop('y', "P124A (V)")
    x = kwargs.pop('x', "B Field (T)")
    X = df[x]
    Y = df[y]
    # Create plot
    fig,ax = plt.subplots()
    # Put data on plot
    ax.scatter(X,Y,label=filename)
    # Place a legend so user knows what they're looking at.
    ax.legend(loc='best')
    # Advise users on how to operate ginput
    print("\n\nClick the suspected peaks (or valleys) of the Oscillation.")
    print("Press ENTER. When all peaks have been denoted.")
    print("Left click to mark, Right click to unmark.")
    # User marks peaks.
    tuples = plt.ginput(-1, timeout=False)
    # Return tuples to parent function for further shaping
    return tuples

def ElectABOPeaks(df,window,filename,identifier,**kwargs):
    analysisdf = pd.DataFrame()
    y = kwargs.pop('y', "P124A (V)")
    x = kwargs.pop('x', "B Field (T)")
    fpath = kwargs.pop("path","")
    selectregion=kwargs.pop("selectregion",True)

    function="ThirdOrder"
    fitname=y+" "+function+" fit"
    fitsub=fitname+" subtraction"
    #df = df.iloc[window[0]:window[1]]
    # Does fit subtraction for whole window
    df,fig = GFF(df,function,filename=filename,window=window)
    if selectregion:
        print("\n\nClick the suspected peaks of the Oscillation.")
        print("Press ENTER. When all peaks have been denoted.")
        print("Left click to mark, Right click to unmark.")
        # User marks peaks.
        tuples = plt.ginput(-1, timeout=False)
        plt.close()
    else:
        return fig,df

    s = StringIO()
    print(window,file=s, end='')
    wd = s.getvalue()

    td = Extract(tuples)
    analysisdf[wd] = td
    return analysisdf



def ingestData(file='stats.csv'):
    with open(file, 'r') as f:
        stats = pd.read_csv(f)
    return stats


def main():
    # Fetch the stats data - it'll be useful for later
    stats = ingestData()
    # Grab measurement numbers from stats.csv
    nums = stats["Measurement"]
    # Look if user has analyzed any previous data.
    try:
        # If previous data exists, then read it in.
        with open('res.csv','r') as f:
            results= pd.read_csv(f)

        cols = results.columns.to_list()
    except FileNotFoundError:
        # If previous data doesn't exist, then create data-structure.
        results = pd.DataFrame()
        cols = []

    for NUM in nums:
        """
            This for loop is what enables the user to elect
            peaks and then save their elected peak's x-coordinates
            to a csv file
        """
        if str(NUM) in cols:
            # If we've already analyzed this data, then don't ask the user to re-select it.
            # INSTEAD: THE USER MUST DELETE THE COLUMN IN THE "res.csv" THAT THEY DISAPPROVE OF.
            continue

        # The datafiles are named "squ-3-dev-31Aug23-NUM.dat
        #   So, in a for loop, I create a file-name-scheme variable that
        #   updates after I've selected the peaks I'm interested in
        FILENAMESCHEME=f"tristan_ringsgate-{NUM}.dat"
        
        # ElectPeaks returns a list of tuples like [(x1, y1), ... ,(xn, yn)],
        tuples = ElectPeaks(FILENAMESCHEME)
        # Unpack the tuples in a syntactically dense manner.
        X,Y=[[x for x,y in tuples],
            [y for x,y in tuples]]
        # Add an additional column into the data structure that are the x-coordinates of the valleys
        results = pd.concat([results,pd.DataFrame({NUM:X})],axis=1)

        # Write that data structure to a file after each loop execution, so the most a user will loose upon an error
        #   Is a single measurement's worth of peak-election.
        with open('res.csv','w') as f:
            results.to_csv(f, index=False)

    # Return True per PEP standards.
    return True
if __name__ == '__main__':
    main()
