import pandas as pd
import numpy as np
import pylab as pl
import ASCII_GUI as gui
import matplotlib.pyplot as plt
import matplotlib,multiprocessing,os
import peak_elector, utilities, windower, record_handler, ABO, hall, quantum_hall
from scipy.optimize import curve_fit as fit
pd.options.mode.chained_assignment = None 
from io import StringIO
from multiprocessing import Pool
#font = {'size'   : 8,
        #"family":'serif'}
#matplotlib.rc('font', **font)
#plt.rcParams.update({
        #"text.usetex": True})
pl.rcParams['figure.figsize']  = 8.5, 11*.75
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



class MeasurementAnalyzer():
	def __init__(self,meas_rec=None, dat_name_template='tristan_ringsgate-{0}.dat', record=None,y="P124A (V)", x="B Field (T)"):
		self.electRecord(record)
		self.displayTypes()
		print(f"Unique Measurment types in: {self.record_path}: "+', '.join(self.mtypes))
		gui.Announcement("Measurement Breakdown")
		print(self.typesdf)
		self.dat_name_template=dat_name_template
		self.map = {}
		self.x=x
		self.y=y
		

		self.mainloop()

	def mainloop(self):
		d = {
		     "ABO":["Aharonov-Bohm data",self.aboTreatment],
		     "Hall":["Hall Data",self.hallTreatment],
		     "QHall":["Quantum Hall",self.quantumHallTreatment],
		     "ABO Plot":["Plot Aharonov-Bohm data",self.aboPlot],
		     "QHall Plot":["Plot Quantum Hall",self.quantumHallPlot]
		    }
		while True:
			gui.Announcement("Make a selection")
			f = d[gui.dict_selector(d)][1]
			f()
	
	def triageTypes(self,caller):
		# Get the user to map a measurement type to a style of 
		d = dict(zip([i for i in range(len(self.mtypes))],self.mtypes))
		gui.Announcement("Which Measurement Type Will you Analyze?")
		self.triage = gui.select_simple(d)

	def aboTreatment(self):
		r = self._treatmentPreamble("ABO")
		ms = r["Measurement"].values.tolist()
		to_analyze = gui.complex_selector(dict(zip(ms,ms)))
		printable = [str(i) for i in to_analyze]
		print("You will now be analyzing Measurements: "+", ".join(printable))
		f = [self.numToFilename(i) for i in to_analyze]
		cf = {f[i]:r[r["Measurement"]==v]["Correction"].values.tolist()[0]\
			 for i,v in enumerate(to_analyze)}
		print(cf)
		windowwidth = self.getABOWidth()
		ABO.ABOMain(f,cf,subwindowWidth=windowwidth)
	
	def getABOWidth(self):
		print("Input window-width in Tesla")
		suffixs = ["T", "mT", "uT", 'G', "mG"]
		scale = [1,10**-3, 10**-6, 10**-4, 10**-7]
		w = input("Enter window-width Ex: {NUM} {T, mT, uT, G, mG}: ")
		wlisted = w.split(' ')
		if len(wlisted) == 1:
			return wlisted[0]
		elif len(wlisted) == 2:
			num, suffix = wlisted
		else:
			print("Improper input.")
			return self.getABOWidth()
		num = float(num)
		for i,s in enumerate(suffixs):
			if s.lower() == suffix.lower():
				return num*scale[i]

	def numToFilename(self,num):
		return self.dat_name_template.format(num)
		
	def aboPlot(self):		
		pass
	
	def quantumHallPlot(self):
		pass
	def _treatmentPreamble(self,name):
		self.triageTypes(name)
		self.map[name] = self.triage
		r = self.record.record
		r = r[r['Type'] == self.triage]
		return r
		
	def hallTreatment(self):
		r = self._treatmentPreamble("Hall")
		ms = r["Measurement"].values.tolist()
		to_analyze = gui.complex_selector(dict(zip(ms,ms)))
		f = [self.numToFilename(i) for i in to_analyze]
		I = r["I_ac"]
		for i,v in enumerate(f):
			n = int(to_analyze[i])
			hall.hall(v, r[r["Measurement"]==n]["I_ac"].values.tolist()[0],selx=True)
		

	def quantumHallTreatment(self):
		e = 1.60217663*10**-19
		h = 6.62607015*10**-34
		Rk = h/(e**2)
		RName = "Klitzing Resistance"
		r = self._treatmentPreamble("QHall")
		ms = r["Measurement"].values.tolist()
		to_analyze = gui.complex_selector(dict(zip(ms,ms)))
		f = [self.numToFilename(i) for i in to_analyze]
		for i,fn in enumerate(f):
			n = int(to_analyze[i])
			m = r[r["Measurement"]==n]
			I = m["I_ac"].values.tolist()[0]
			c = m["Correction"].values.tolist()[0]
			df = utilities.ReadDatFile(fn, correction=c, y=self.y)
			df[RName] = (df[self.y].values/I)/Rk
			dpts = []
			addition = [None,None]
			fig,ax = None,None
			while len(addition) ==2:
				addition = windower.selectWindower(df,self.x,\
						self.y,clickpoints=True, fig=fig,\
						ax=ax)
				dpts += addition
				fig,ax = plt.subplots(2)
				for x in dpts:
					ax[0].axvline(x=x, color='red', linestyle='--')
				
				
					
				
				
			
			
		
		pass
	
	def electRecord(self, record=None):
		# asciigui interface: allow user to select measurement record.
		# Record is a RecordHandler class defined in record_handler.py
		if record is None:
			potential_records = sniff_dirs('record', exclude=[])
			choices = [i for i in range(len(potential_records))]
			gui.Announcement("Select Record")
			self.record_path = gui.select_simple(dict(\
				    zip(choices,potential_records)))
		else:
			self.record_path = record
		self.record = record_handler.Record(self.record_path)
		
	
	def displayTypes(self):
		m = "Measurement"
		self.mtypes = self.record.get_types()
		# I want an n-column of types
		# an m-row of measurments
		r = self.record.record
		tempdf = pd.DataFrame()
		for t in self.mtypes:
			tempdf = pd.concat([tempdf, \
				pd.DataFrame({t:r[r["Type"]==t][m].values.tolist()})]\
				,ignore_index=True,axis=1)
		c = tempdf.columns.to_list()
		self.typesdf = tempdf.rename(columns=dict(zip(c,self.mtypes)))

	def read_measurement_record(self):
		with open(self.record_path, 'r') as f:
			self.record = pd.read_csv(f)


def sniff_dirs(search:str, exclude=[]):
	# returns a list of stuff that contains search
	items = sorted(os.listdir('.'))
	match = []
	for i in items:
		if search.lower() in i.lower():
			match.append(i)
	return match
	
	

def main():
	MA = MeasurementAnalyzer(record='Measurement_Record.csv')

if __name__ == '__main__':
	main()

