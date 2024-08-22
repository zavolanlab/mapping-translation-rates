#!/usr/bin/env python
'''
script uses binned depletion curve files
to determine the global elongation speed
via the SL method
as established by Ingolia et. al. 2011
'''
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


files = ['output_HepG2_normal.txt','output_HepG2_HS.txt','output_Huh7_normal.txt']



#parameters
omission = 120
SL_binsize = 15
threshold = 0.5


'''
linear fit function
at point x with slope m*m and intercept b
CAUTION: m>=0 enforced!
'''
def linear(x, m, b):
	return m*m*x + b



'''
BEGIN SCRIPT
'''
fig, ax = plt.subplots()

plt.xlabel('time difference [s]')
plt.ylabel('starting location [nt]')

plt.xlim([70,200])
plt.ylim([100,2000])

for filename in files:
	df = pd.read_csv(filename,sep="\t",header=None,low_memory=False)

	df = df.transpose()
	df.columns = df.iloc[0]
	df = df[1:]
	df = df.drop(columns=['t'])

	reference_profile = df['0'].to_numpy()

	columns = df.columns.tolist()
	columns.remove('0')

	locations = {}
	for t in columns:
		profile = df[t].to_numpy()

		profile = profile / reference_profile

		for i in range(omission//SL_binsize,len(profile)-1):
			if profile[i] < threshold and profile[i+1] >= threshold:
				locations[t] = (i + 0.5) * SL_binsize
				break

	x_data = np.array([*locations.keys()],dtype=float)
	y_data = np.array([*locations.values()],dtype=float)


	if len(y_data) > 1:
		popt, _ = curve_fit(linear, x_data, y_data, maxfev=2000)

		condition = filename.replace("output_", "").replace(".txt", "")

		plt.scatter(x_data,y_data,label=condition.replace("_"," "))

		function_values = linear(x_data,*popt)
		plt.plot(x_data,function_values,label="linear fit, slope = %.2f nt/s"%(popt[0]*popt[0]))

		plt.legend(loc="upper left")

	plt.savefig('fits.pdf', bbox_inches='tight')

