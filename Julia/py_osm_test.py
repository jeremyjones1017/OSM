import numpy as np
import csv
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from astropy.io import ascii

def read_filters(use_filts,filt_dir,cwl,wav):
	"""Reads and stores the data from the filter curves to be used.
	Inputs:
	use_filts
		A list of the names of the filters used by the photometry
	filt_dir
		The directory in which the filter curves are stored
	cwl
		A dictionary with the values for the central wavelengths
	wav
		An array of wavelength values used by the model spectra
	Outputs:
	filt_dict
		A dictionary with the filter transmission curve expressed where 
		the corresponding wavelength array is that being used for the SED
	"""
	filt_dict=dict()
	for i in range(len(use_filts)):
		#print '========================================'
		#print use_filts[i]
		#print '========================================'
		this_wav=[]
		this_tran=[]
		with open(filt_dir+use_filts[i]+'.txt') as input:	#This opens and records the wavelength and transmission vectors
			try:
				input_reader=csv.reader(input,delimiter='\t')	#If the data are tab-delimited
				input_reader.next()
				for line in input_reader:
					this_wav.append(float(line[0])/1e8)
					this_tran.append(float(line[1]))
			except:
				input_reader=csv.reader(input,delimiter=' ')	#If the data are space-delimited
				input_reader.next()
				for line in input_reader:
					this_wav.append(float(line[0])/1e8)
					this_tran.append(float(line[1]))
		#plt.plot(this_wav,this_tran,'b')
		int_trans=np.zeros(len(wav))
		f=interp1d(this_wav,this_tran)	#This interpolates the saved filter curve to the appropriate wavelength scale
		for j in range(len(wav)):
			if wav[j] >= this_wav[0] and wav[j] <= this_wav[-1]:
				int_trans[j]=f(wav[j])
		#plt.plot(wav,int_trans,'r--')
		#plt.show()
		abs_wav_minus_cwl=abs(wav - cwl[use_filts[i]])
		cwl_trans=int_trans[np.where(abs_wav_minus_cwl == min(abs_wav_minus_cwl))]
		int_trans=int_trans/cwl_trans[0]
		filt_dict[use_filts[i]]=int_trans
	return filt_dict
def read_cwlzpf(cwlzpf_file):
	data=ascii.read(cwlzpf_file)
	filt=list(data['waveband'])
	cwl=list(data['central_wavelength_in_cm'])
	zpf=list(data['zero_point_flux_in_erg/s/cm^2/cm'])
	cwl_dict=dict()
	zpf_dict=dict()
	for i in range(len(filt)):
		cwl_dict[filt[i]]=float(cwl[i])
		zpf_dict[filt[i]]=float(zpf[i])
	return cwl_dict,zpf_dict


def read_phot(phot_inp):
	phot_data=dict()
	use_filts=[]
	with open(phot_inp,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			if line[0][0] != '#':
				phot_data[line[0]]=[line[1],line[2]]
				use_filts.append(line[0])
	return phot_data,use_filts

phot_inp="C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/HD192640/HD192640.phot"
phot_data,use_filts=read_phot(phot_inp)

for i in phot_data:
	print i,phot_data[i]

'''
cwl,zpf=read_cwlzpf("C:/Users/Jeremy/Dropbox/Programing/Astars/Band_Passes/cwlzpf.txt")
#use_filts=["2massH","2massJ","2massK","cousinsI","cousinsR","F1565","F1965","F2365","F2740","johnsonB","johnsonH","johnsonI","johnsonJ","johnsonK","johnsonR","johnsonR","johnsonU","johnsonV","stromgrenb","stromgrenu","stromgrenv","stromgreny","wes15N","wes15W","wes18","wes22","wes25","wes33"]
use_filts=[]
for i in cwl:
	use_filts.append(i)
phx_wav=(np.arange(25500)+500.)*1e-8 #This defines the wavelength array
filt_dict=read_filters(use_filts,"C:/Users/Jeremy/Dropbox/Programing/Astars/Band_Passes/",cwl,phx_wav)
for i in filt_dict:
	print i,filt_dict[i],sum(filt_dict[i])
'''


