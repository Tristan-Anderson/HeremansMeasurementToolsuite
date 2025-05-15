from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

p = r"C:\Users\Tristan\Desktop\work\Research\heremans\MonarkMeasurements\20250203Cooldown1AmpOhmn\Mixer"
filesn = [17,18,19,20,21]
bend = [19,17]
sym = [20,21,18]
key = {19:"10 um Bend", 20:"10 um Symmetric (+)",21:"10 um Symmetric (-)",17:"100 um Bend", 18:"100 um Symmetric"}
files = [f"{p}\\{i}.csv" for i in filesn]
dfs = []
for i in files:
    with open(i, 'r') as f:
        dfs.append(pd.read_csv(f))
v = "lockin_1 - X"
for d in range(len(dfs)):
    df = dfs[d]
    df["I_dc"] = df[v]*100*10**-6
    #df["Rx"] = df["lockin_1 - X"]/(100*10**-9) # All taken at 100nA
    dfs[d] = df

pfiles = [f"{p}\\{i}.png" for i in filesn]
for i, fn in enumerate(pfiles):
    df = dfs[i]
    x = df["I_dc"]
    y = df[v]
    plt.plot(x,y,label=key[filesn[i]])
plt.legend(loc='Best')
plt.savefig(f"{p}\\test.png")



