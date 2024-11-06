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

