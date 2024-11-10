import pandas as pd
import numpy as np
import pylab as pl
import ASCII_GUI as gui
import matplotlib.pyplot as plt
import matplotlib,multiprocessing,os
import peak_elector, utilities, windower, record_handler
from scipy.optimize import curve_fit as fit
pd.options.mode.chained_assignment = None 
from io import StringIO
from multiprocessing import Pool
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

class MeasurementAnalyzer():
	def __init__(self,meas_rec=None, dat_name_template='tristan_ringsgate-{NUM}.dat', record=None):
		self.electRecord(record)
		self.displayTypes()
		print(f"Unique Measurment types in: {self.record_path}: "+', '.join(self.mtypes))
		gui.Announcement("Measurement Breakdown")
		print(self.typesdf)
		self.dat_name_template=dat_name_template
		self.mainloop()

	def mainloop()
	
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
			tempdf = pd.concat([tempdf, pd.DataFrame({t:r[r["Type"]==t][m].values.tolist()})],ignore_index=True,axis=1)
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

