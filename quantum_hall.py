import pandas as pd
import numpy as np
import pylab as pl
import utilities, windower
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit

### General function to extract carrier densities from measurements
### Selects single region
def quantum_hall(
		fn, # Df containing data
		r,  # Df (record)
		n,  # Measurement number
		x='B Field (T)', 
		y='P124A (V)', 
		f=None):
	"""
	    Function which handles Hall Data. Requests singular file, with paramaters: x,y, and selx.
	    User can elect x-selection capabilities on call()
	"""
	e = 1.60217663*10**-19
	h = 6.62607015*10**-34
	Rk = h/(e**2)
	RName = "Klitzing Resistance"
	m = r[r["Measurement"]==n]
	I = m["I_ac"].values.tolist()[0]
	c = m["Correction"].values.tolist()[0]
	df = utilities.ReadDatFile(fn, correction=c, y=y)
	df[RName] = (df[y].values/I)/Rk
	dpts = []
	addition = [None,None]
	fig,ax = None,None
	while len(addition) == 2:
		addition = windower.selectWindower(df,x,\
				RName,clickpoints=True, fig=fig,\
				ax=ax,filename=fn)
		dpts += addition
	else:
		if len(addition) == 0:
			pass
		else:
			cols = [i for i in range(len(dpts))]
			results = {i:dpts[i] for i in cols}
			rdf = pd.DataFrame(results)
			print(rdf)
			with open("QH_Record_"+fn.split('.dat')[0]+'.csv','w') as recordfile:
				rdf.to_csv(recordfile)		
			plt.close()
		return dpts

def quantum_hall_density(
		df,**kwargs):
	# elect peaks
	#fwindow = quantum_hall(fn,r,n,x,y,f)
	fwindow = [[np.float64(5.919787309021314), np.float64(5.951316240043938)],
				[np.float64(4.05402684335543), np.float64(4.320053190456367)],
				[np.float64(2.4518113425027765), np.float64(2.5620817679257355)],
				[np.float64(1.9986063638096052), np.float64(2.1327172798093783)],
				[np.float64(1.7514540002187564), np.float64(1.8376500836898595)],
				[np.float64(1.508868647336885), np.float64(1.6024654046543128)],
				[np.float64(1.3672556955612944), np.float64(1.4220942647970858)],
				[np.float64(1.2204778776046354), np.float64(1.279646679544388)],
				[np.float64(1.1261801940173943), np.float64(1.1567346220107373)],
				[np.float64(0.9525360758294584), np.float64(0.9757550893931559)],
				[np.float64(0.8832747348394816), np.float64(0.9063223111212964)],
				[np.float64(0.7340236706265857), np.float64(0.744162935078455)],
				[np.float64(0.6588254848500176), np.float64(0.6665119511714761)],
				[np.float64(0.4808505908992973), np.float64(0.4865753650369775)],
				[np.float64(-0.6995553068472785), np.float64(-0.6850621376063435)],
				[np.float64(-0.7896086407008956), np.float64(-0.7726937406626075)],
				[np.float64(-0.8447447638869379), np.float64(-0.8273642872986182)],
				[np.float64(-0.8436515207536843), np.float64(-0.827760254826997)],
				[np.float64(-0.9070521513446363), np.float64(-0.8841520189883243)],
				[np.float64(-0.9799382897937754), np.float64(-0.9558224077757485)],
				[np.float64(-1.059540680509004), np.float64(-1.0228031053715747)],
				[np.float64(-1.1593416001804187), np.float64(-1.1262969465620634)],
				[np.float64(-1.2771501268106564), np.float64(-1.2242193380705555)],
				[np.float64(-1.4177674854518971), np.float64(-1.3689016888369652)],
				[np.float64(-1.6027130160232477), np.float64(-1.5193830199421774)],
				[np.float64(-1.835977648789938), np.float64(-1.7522853483780805)],
				[np.float64(-2.126300360097899), np.float64(-1.9934533566004393)],
				[np.float64(-2.5699061391001234), np.float64(-2.455688920269014)],
				[np.float64(-3.2223861717230724), np.float64(-2.9479507919196104)],
				[np.float64(-4.322757065904571), np.float64(-4.060136399094949)],
				[np.float64(-5.951689668680859), np.float64(-5.925301502539693)]]
	fwindowp = []
	for l in fwindow:
		fwindowp += l
	bnuwindow = [max(i)-min(i) for i in fwindow]
	print(fwindow)
	print(bnuwindow)
	plt.close('all')

	# Get Resistance
	x = kwargs.pop('x',"B Field (T)")
	y = kwargs.pop('y',"PV124A (V)")
	R = kwargs.pop('R',"Resistance")
	I = kwargs.pop("I",0)
	e = 1.60217663 * 10 ** -19
	h = 6.62607015 * 10 ** -34
	Rk = h / (e ** 2)
	RName = "Klitzing Resistance"
	df[R] = (df[y].values / I)
	df[RName] = df[R]/ Rk

	# Make figure
	fig,ax = plt.subplots(2)

	ax[0].plot(df[x], df[R].values/1000, color='blue', label="Resistance")
	ax[0].set_ylabel(R + " (kOhms)")
	ax[1].plot(df[x], df[RName], color='blue', label="1/v")
	ax[1].set_ylabel(RName + " (Normalized)")

	#
	figxrange = max(df[x])-min(df[x])
	figxstep = figxrange/20
	Rfigrange = (max(df[R])-min(df[R]))/1000
	Vfigstep = Rfigrange/20

	Rfigrange2 = max(df[RName])-min(df[RName])
	Vfigstep2 = Rfigrange2/20

	# For every window that was elected
	for idx,i in enumerate(fwindow):
		# Unpack window and get important values.
		even = idx%2
		if even:
			xloc = max(df[x])-5*figxstep
		else:
			xloc = min(df[x])+5*figxstep
		print(i)
		x1,x2 = i
		center = bnuwindow[idx] 	# Is this the Bnu value?
		cut = df[(df[x]>x1)&(df[x]<x2)]
		# Get the Nu value: mean nu over the cut window
		v_value = np.mean(1/cut[RName].values)
		print(center)
		print(v_value)
		#Ns = center*v_value*2.418*10**10
		# eB/h*v  = Ns
		Ns = (e*center/h)*v_value

		# Get the Resistance in kOhms
		yvals = cut[R].values/1000
		ymean = np.mean(yvals)
		# Plot the kOhm markers.
		ax[0].axhline(y=ymean, color='red', linestyle='--')

		ax[0].text(xloc, ymean+Vfigstep/2, f"v: {v_value:0.1f} N_s: {Ns:0.2E}")

		# Get 1/v values.
		yvals = cut[RName].values
		ymean = np.mean(yvals)
		# Plot the 1/v values.
		ax[1].axhline(y=ymean, color='red', linestyle='--')
		ax[1].text(xloc, ymean, f"v: {v_value:0.1f} N_s: {Ns:0.2E}")

	plt.show()

def qh_agile(
		df, # Df containing data
		**kwargs):
	x = kwargs.pop('x',"B Field (T)")
	y = kwargs.pop('y',"PV124A (V)")
	fn = kwargs.pop('fn',"")
	I = kwargs.pop("I",0)
	e = 1.60217663 * 10 ** -19
	h = 6.62607015 * 10 ** -34
	Rk = h / (e ** 2)
	RName = "Klitzing Resistance"
	df[RName] = (df[y].values / I) / Rk
	dpts = []
	addition = [None, None]
	fig, ax = None, None
	try:
		while len(addition) == 2:
			addition = windower.selectWindower(df, x, \
										   RName, clickpoints=True, fig=fig, \
										   ax=ax, filename=fn,dpts=dpts)
			print(addition)
			dpts.append(addition)
		else:
			if len(addition) == 0:
				pass
			else:
				cols = [i for i in range(len(dpts))]
				results = {i: dpts[i] for i in cols}
				rdf = pd.DataFrame(results)
				print(rdf)
				with open("QH_Record_" + fn.split('.dat')[0] + '.csv', 'w') as recordfile:
					rdf.to_csv(recordfile)
				plt.close()
			return dpts
	except KeyboardInterrupt:
		cols = [i for i in range(len(dpts))]
		results = {i: dpts[i] for i in cols}
		rdf = pd.DataFrame(results)
		with open("QH_Record_" + fn.split('.dat')[0] + '.csv', 'w') as recordfile:
			rdf.to_csv(recordfile)


def quantum_hall_density2(df, **kwargs):
	# elect peaks
	# fwindow = quantum_hall(fn,r,n,x,y,f)
	fwindow = [[np.float64(5.92251075660082), np.float64(5.9512951508306156)],
[np.float64(4.050777224424044), np.float64(4.324873292668928)],
[np.float64(2.932421910357532), np.float64(3.2263224879408297)],
[np.float64(2.44929450480496), np.float64(2.56837202811305)],
[np.float64(1.991510611993114), np.float64(2.1394619204849947)],
[np.float64(1.7485500181879283), np.float64(1.8413338956443897)],
[np.float64(1.5032026799245328), np.float64(1.6087740386865197)],
[np.float64(1.362627184146566), np.float64(1.4251750685531237)],
[np.float64(1.2151829290594724), np.float64(1.2824332054819272)],
[np.float64(1.123230063548885), np.float64(1.16271310128682)],
[np.float64(1.0205538894505448), np.float64(1.063401058152582)],
[np.float64(0.9525057321650574), np.float64(0.9815266983146991)],
[np.float64(0.8779867085230197), np.float64(0.9116605193129199)],
[np.float64(0.8287286117054491), np.float64(0.8495729061151908)],
[np.float64(0.7694443020252129), np.float64(0.7948273681120008)],
[np.float64(0.7321627258610155), np.float64(0.7479060875494532)],
[np.float64(0.6870722775864424), np.float64(0.7052237840397715)],
[np.float64(0.6558067056540491), np.float64(0.6672741232425299)],
[np.float64(0.6201002852841324), np.float64(0.6320039720688786)],
[np.float64(0.566887696300791), np.float64(0.5783784805965811)],
[np.float64(0.5400958344676773), np.float64(0.5511402680972832)],
[np.float64(0.590378104797533), np.float64(0.6043330173165112)],
[np.float64(0.4768986463029502), np.float64(0.48951260638010397)],
[np.float64(0.5197994470005134), np.float64(0.5293147666286846)],
[np.float64(0.4451989778846289), np.float64(0.45406884899630073)],
[np.float64(0.38831576422701813), np.float64(0.3979169430773323)],
[np.float64(-5.951964520703563), np.float64(-5.921676005062481)],
[np.float64(-4.32807400436052), np.float64(-4.054221281298801)],
[np.float64(-3.227131418832552), np.float64(-2.9467116148274397)],
[np.float64(-2.5707366181931874), np.float64(-2.4565758489213576)],
[np.float64(-2.143933599213037), np.float64(-1.9891681615409837)],
[np.float64(-1.8392366185744677), np.float64(-1.7491202486712365)],
[np.float64(-1.6071644333330752), np.float64(-1.509981546046771)],
[np.float64(-1.4264147120949848), np.float64(-1.365085750779719)],
[np.float64(-1.2814838269587543), np.float64(-1.218751279668564)],
[np.float64(-1.1573178738815508), np.float64(-1.1242598423519385)],
[np.float64(-1.064009629314018), np.float64(-1.0254231103536429)],
[np.float64(-0.9783490446057207), np.float64(-0.9526305710678572)],
[np.float64(-0.9098150230695901), np.float64(-0.876180869399428)],
[np.float64(-0.8439239990084976), np.float64(-0.8274218773480406)],
[np.float64(-0.7927804809561804), np.float64(-0.7724738897882989)],
[np.float64(-0.7461578824030676), np.float64(-0.7283919514273539)],
[np.float64(-0.7036232204965088), np.float64(-0.6832277546833587)],
[np.float64(-0.6663520275348058), np.float64(-0.6538230059733638)],
[np.float64(-0.6318931477569399), np.float64(-0.6207467949924232)],
[np.float64(-0.5981776015393211), np.float64(-0.5928954823757171)],
[np.float64(-0.5737173596443739), np.float64(-0.5626894631297227)],
[np.float64(-0.5250796761823604), np.float64(-0.5139212611844908)],
[np.float64(-0.4817535129284899), np.float64(-0.47598164941556137)],
[np.float64(-0.4500174357971565), np.float64(-0.4423915070479659)],
[np.float64(-0.3900786406391436), np.float64(-0.3858943720728875)],
[np.float64(-0.34846175181217875), np.float64(-0.34291649889419706)]]

	dump = kwargs.pop('dump', '')
	fwindowp = []
	for l in fwindow:
		fwindowp += l
	plt.close('all')

	# Get Resistance
	x = kwargs.pop('x', "B Field (T)")
	y = kwargs.pop('y', "PV124A (V)")
	R = kwargs.pop('R', "Resistance")
	I = kwargs.pop("I", 0)
	e = 1.60217663 * 10 ** -19
	h = 6.62607015 * 10 ** -34
	Rk = h / (e ** 2)
	RName = "Klitzing Resistance"
	df[R] = (df[y].values / I)
	df[RName] = df[R] / Rk
	# Make figure

	# For every window that was elected
	for idx, i in enumerate(fwindow):
		x1,x2 = min(i), max(i)
		if x1 <0:
			x2,x1 = min(i), max(i)
		if x1 >0:
			cut = df[(df[x]>=x1*.95) & (df[x]<=x2*1.05)]
		else:
			cut = df[(df[x] <= x1 * .95) & (df[x] >= x2 * 1.05)]
		fig, ax = plt.subplots(2)

		ax[0].plot(df[x], df[R].values / 1000, color='blue', label="Resistance")
		ax[0].plot(cut[x], cut[R]/1000, color='red', label="Selected Plateau")
		ax[0].set_ylabel(R + " (kOhms)")
		ax[1].plot(cut[x], 1/cut[RName], color='blue', label="v")
		ax[1].set_ylabel(f"Resistance/{RName}")

		center = (x1+x2)/2  # center between values.
		if x1 >0:
			cut = df[(df[x] >= x1) & (df[x] <= x2)]
		else:
			cut = df[(df[x] <= x1) & (df[x] >= x2)]
		ax[1].plot(cut[x], 1/cut[RName], color='red', label="Plateau")
		# Get the Nu value: mean nu over the cut window
		v_value = np.mean(1 / cut[RName].values)
		print("Center: ", center)
		print("Nu Value: ",v_value)
		# Ns = center*v_value*2.418*10**10
		# eB/h*v  = Ns
		Ns = (e * center / h) * v_value

		ax[1].axhline(y=v_value, color='red', linestyle='--',label=f"Nu = {v_value:0.02f}")
		ax[1].axvline(x=center,color='blue',linestyle = '--',label=f"Center x:{center:0.02E}")

		fig.suptitle(f"v: {v_value:0.2f} N_s: {Ns:0.03E}")
		for v in ax:
			v.legend(loc='best',frameon=False)

		# Get 1/v values.
		plt.savefig(f"{dump}{idx}")

def rxx(df,**kwargs):
	# elect peaks
	# fwindow = quantum_hall(fn,r,n,x,y,f)
	fwindow = [[np.float64(5.92251075660082), np.float64(5.9512951508306156)],
			   [np.float64(4.050777224424044), np.float64(4.324873292668928)],
			   [np.float64(2.932421910357532), np.float64(3.2263224879408297)],
			   [np.float64(2.44929450480496), np.float64(2.56837202811305)],
			   [np.float64(1.991510611993114), np.float64(2.1394619204849947)],
			   [np.float64(1.7485500181879283), np.float64(1.8413338956443897)],
			   [np.float64(1.5032026799245328), np.float64(1.6087740386865197)],
			   [np.float64(1.362627184146566), np.float64(1.4251750685531237)],
			   [np.float64(1.2151829290594724), np.float64(1.2824332054819272)],
			   [np.float64(1.123230063548885), np.float64(1.16271310128682)],
			   [np.float64(1.0205538894505448), np.float64(1.063401058152582)],
			   [np.float64(0.9525057321650574), np.float64(0.9815266983146991)],
			   [np.float64(0.8779867085230197), np.float64(0.9116605193129199)],
			   [np.float64(0.8287286117054491), np.float64(0.8495729061151908)],
			   [np.float64(0.7694443020252129), np.float64(0.7948273681120008)],
			   [np.float64(0.7321627258610155), np.float64(0.7479060875494532)],
			   [np.float64(0.6870722775864424), np.float64(0.7052237840397715)],
			   [np.float64(0.6558067056540491), np.float64(0.6672741232425299)],
			   [np.float64(0.6201002852841324), np.float64(0.6320039720688786)],
			   [np.float64(0.566887696300791), np.float64(0.5783784805965811)],
			   [np.float64(0.5400958344676773), np.float64(0.5511402680972832)],
			   [np.float64(0.590378104797533), np.float64(0.6043330173165112)],
			   [np.float64(0.4768986463029502), np.float64(0.48951260638010397)],
			   [np.float64(0.5197994470005134), np.float64(0.5293147666286846)],
			   [np.float64(0.4451989778846289), np.float64(0.45406884899630073)],
			   [np.float64(0.38831576422701813), np.float64(0.3979169430773323)],
			   [np.float64(-5.951964520703563), np.float64(-5.921676005062481)],
			   [np.float64(-4.32807400436052), np.float64(-4.054221281298801)],
			   [np.float64(-3.227131418832552), np.float64(-2.9467116148274397)],
			   [np.float64(-2.5707366181931874), np.float64(-2.4565758489213576)],
			   [np.float64(-2.143933599213037), np.float64(-1.9891681615409837)],
			   [np.float64(-1.8392366185744677), np.float64(-1.7491202486712365)],
			   [np.float64(-1.6071644333330752), np.float64(-1.509981546046771)],
			   [np.float64(-1.4264147120949848), np.float64(-1.365085750779719)],
			   [np.float64(-1.2814838269587543), np.float64(-1.218751279668564)],
			   [np.float64(-1.1573178738815508), np.float64(-1.1242598423519385)],
			   [np.float64(-1.064009629314018), np.float64(-1.0254231103536429)],
			   [np.float64(-0.9783490446057207), np.float64(-0.9526305710678572)],
			   [np.float64(-0.9098150230695901), np.float64(-0.876180869399428)],
			   [np.float64(-0.8439239990084976), np.float64(-0.8274218773480406)],
			   [np.float64(-0.7927804809561804), np.float64(-0.7724738897882989)],
			   [np.float64(-0.7461578824030676), np.float64(-0.7283919514273539)],
			   [np.float64(-0.7036232204965088), np.float64(-0.6832277546833587)],
			   [np.float64(-0.6663520275348058), np.float64(-0.6538230059733638)],
			   [np.float64(-0.6318931477569399), np.float64(-0.6207467949924232)],
			   [np.float64(-0.5981776015393211), np.float64(-0.5928954823757171)],
			   [np.float64(-0.5737173596443739), np.float64(-0.5626894631297227)],
			   [np.float64(-0.5250796761823604), np.float64(-0.5139212611844908)],
			   [np.float64(-0.4817535129284899), np.float64(-0.47598164941556137)],
			   [np.float64(-0.4500174357971565), np.float64(-0.4423915070479659)],
			   [np.float64(-0.3900786406391436), np.float64(-0.3858943720728875)],
			   [np.float64(-0.34846175181217875), np.float64(-0.34291649889419706)]]

	dump = kwargs.pop('dump', '')
	fwindowp = []
	for l in fwindow:
		fwindowp += l
	plt.close('all')

	# Get Resistance
	x = kwargs.pop('x', "B Field (T)")
	y = kwargs.pop('y', "PV124A (V)")
	yxy = kwargs.pop('yxy','lockin_1 - X')
	R = kwargs.pop('R', "Resistance")
	I = kwargs.pop("I", 0)
	e = 1.60217663 * 10 ** -19
	h = 6.62607015 * 10 ** -34
	Rk = h / (e ** 2)
	RName = "Klitzing Resistance"
	df[R] = ((df[y].values)/ I)
	df[R+'xy'] = ((df[yxy].values*-1)/I)
	df[RName] = df[R+'xy'].values / Rk
	# Make figure

	# For every window that was elected
	for idx, i in enumerate(fwindow):
		x1, x2 = min(i), max(i)
		if x1 < 0:
			x2, x1 = min(i), max(i)
		if x1 > 0:
			cut = df[(df[x] >= x1 * .95) & (df[x] <= x2 * 1.05)]
		else:
			cut = df[(df[x] <= x1 * .95) & (df[x] >= x2 * 1.05)]
		fig, ax = plt.subplots(2)

		ax[0].plot(df[x], df[R].values / 1000, color='blue', label="Resistance")
		ax[0].plot(cut[x], cut[R] / 1000, color='red', label="Selected Plateau")
		ax[0].set_ylabel(R + " (kOhms)")
		ax[1].plot(cut[x], cut[R], color='blue', label="Resistance")
		ax[1].set_ylabel(R)

		center = (x1 + x2) / 2  # center between values.
		if x1 > 0:
			cut = df[(df[x] >= x1) & (df[x] <= x2)]
		else:
			cut = df[(df[x] <= x1) & (df[x] >= x2)]
		ax[1].plot(cut[x], cut[R], color='red', label="Plateau")
		# Get the Nu value: mean nu over the cut window
		v_value = np.mean(1/cut[RName].values)
		print("Center: ", center)
		Ns = center*v_value*2.418*10**10
		# eB/h*v  = Ns

		ax[1].axvline(x=center, color='blue', linestyle='--', label=f"Center x:{center:0.02E}")

		fig.suptitle(f"v: {v_value:0.03f}")
		for v in ax:
			v.legend(loc='best', frameon=False)

		# Get 1/v values.
		plt.savefig(f"{dump}{idx}")


if __name__ == '__main__':
	#f = [f'tristan_ringsgate-{i}.dat' for i in [17,14]]
	import pylab as pl

	pl.rcParams['figure.figsize'] = 8, 8
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
			df = pd.read_csv(fi)
		df = df.dropna()
		df["B (T)"] = df["AMI_430_Z_Axis_QCoDeS - B"]
		y = "lockin_2 - X"
		df[y] = df[y].values*-1
		quantum_hall_density2(df, selx=True, x="B (T)",
							  y="lockin_1 - X",I=100*10**-9,fn="CombinedHall.csv", dump=path+"\\Rxy")
		rxx(df, selx=True, x="B (T)",
							  y="lockin_2 - X",I=100*10**-9,fn="CombinedHall.csv", dump=path+"\\Rxx")
		#qh_agile(df, x="B (T)", y="lockin_2 - X",selx=True,I=100*10**-9)
