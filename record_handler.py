import numpy as np
import pandas as pd
from utilities import FunctionInputHandler as ReadFile

class Record():

	def __init__(self,path):
		self.record, self.path = ReadFile(path)
		self.cols = self.record.columns.to_list()
		#self.cols = [s.lower() for s in self.cols]
		self.required =['Measurement', 'I_ac', 'Type', 'V+', 'V-', "I+", "I-"]
		self.optional =["Vg", 'I_dc']
		for r in self.required:
			try:
				assert r in self.cols
			except AssertionError:
				raise AssertionError(f"{r}, a required column was not found in the columns")
		measurements = self.record
		
	
	def get_types(self):
		return self.get_unique(self.required[2])

	def get_unique(self,*args):
		# Return unique values (singular, or joint) found in column(s) 
		unique = []
		largs = len(args)
		
		for i,row in self.record.iterrows():
			if largs == 1:
				element = row[args[0]]
			else:
				element = [row[a] for a in args]
			if element not in unique:
				unique.append(element)
		return unique
	
	def get_unique_contacts(self):
		return self.get_unique('V+', 'V-')

	def get_unique_currents(self):
		return self.get_unique('I+', 'I-')
		
