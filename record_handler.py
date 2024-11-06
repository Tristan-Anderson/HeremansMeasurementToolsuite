import numpy as np
import pandas as pd

class Record():

	def __init__(self,path):
		with open(path, 'r') as f:
			self.record = pd.read_csv(f)
		self.path = path
		self.cols = self.record.columns.to_list()
		self.required =['Measurement', 'I_ac', 'Type', 'V+', 'V-', "I+", "I-"]
		self.optional =["Vg", 'I_dc']
		for r in self.required:
			assert r in self.cols

	def get_types(self):
		return self.get_unique(self.required[2])
		
