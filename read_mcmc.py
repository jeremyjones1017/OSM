import numpy as np
from math import *
import csv
from time import *
from scipy.spatial import ConvexHull
from inside import *
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import *
from matplotlib import pyplot as plt
from matplotlib import image as mpimg
import matplotlib.cm as cm
import pyfits
import os
import platform
import multiprocessing
from random import *
import pandas as pd
from pandas.tools.plotting import scatter_matrix
from astropy.io import ascii

def main():
	star='HD31295'
	model='consolidated_vZ'
	mesa_met='_M-0.28'
	star_dir='C:/Users/Jeremy/Dropbox/Python/Astars/Stars/'+star+'/'
	model_dir='C:/Users/Jeremy/Dropbox/Python/Astars/Stars/'+star+'/'+model+'/'
	inp_file=model_dir+star+'.mcmc'
	age_file=model_dir+star+mesa_met+'.ages'
	out_file=model_dir+star+'.results'
	
	print '--------------------------------'
	print star, model
	print '--------------------------------'
	
	run_ages=True
	
	nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff=read(inp_file)
	if run_ages:
		anums,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg=read_ages(age_file)
	else:
		anums=[]
		age=[]
		mass=[]
		omg_init=[]
		L_bol=[]
		L_app=[]
		R_avg=[]
		R_p=[]
		T_e=[]
		T_avg=[]
		lg_p=[]
		lg_e=[]
		lg_avg=[]
	plotting(nums,anums,chi2s,R_e,V_e,inc,T_p,pa,acc,model_dir,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg,run_ages,cutoff,out_file)

def read(inp_file,flag=''):
	nn=0
	R_e=[]
	V_e=[]
	inc=[]
	T_p=[]
	pa=[]
	chi2s=[]
	if flag == '':
		acc=[1.]
		acc_Re=[1.]
		acc_Ve=[1.]
		acc_inc=[1.]
		acc_Tp=[1.]
		acc_pa=[1.]
	if flag == 'ages':
		acc=[]
		acc_Re=[]
		acc_Ve=[]
		acc_inc=[]
		acc_Tp=[]
		acc_pa=[]
	this_Re=[]
	this_Ve=[]
	this_inc=[]
	this_Tp=[]
	this_pa=[]
	nums=[]
	
	with open(inp_file,'r') as input:
		input_reader=csv.reader(input,delimiter='\t')
		input_reader.next()
		for line in input_reader:
			if line == []:
				print 'blank line'
			elif line[0] != 'Num':
				nums.append(int(line[0]))
				R_e.append(float(line[2]))
				V_e.append(float(line[3]))
				inc.append(float(line[4]))
				T_p.append(float(line[5]))
				pa.append(float(line[6]))
				this_Re.append(float(line[2]))
				this_Ve.append(float(line[3]))
				this_inc.append(float(line[4]))
				this_Tp.append(float(line[5]))
				this_pa.append(float(line[6]))
				chi2s.append(float(line[1]))
				if 'consolidated' not in inp_file:
					if int(line[0]) > 0:
						rate_Re_100=float(line[10])
						rate_Ve_100=float(line[11])
						rate_inc_100=float(line[12])
						rate_Tp_100=float(line[13])
						rate_pa_100=float(line[14])
						rate_Re=float(line[15])
						rate_Ve=float(line[16])
						rate_inc=float(line[17])
						rate_Tp=float(line[18])
						rate_pa=float(line[19])
					scale=[float(line[20]),float(line[21]),float(line[22])*pi/180.,float(line[23]),float(line[24])*pi/180.]
					scale=np.array(scale)
					if line[7] == '1.0':
						nn=0
						if int(line[0]) > 0:
							acc.append(1.)
							if line[9] == 'R_e':
								acc_Re.append(1.)
							if line[9] == 'V_e':
								acc_Ve.append(1.)
							if line[9] == 'inc':
								acc_inc.append(1.)
							if line[9] == 'T_p':
								acc_Tp.append(1.)
							if line[9] == 'pa':
								acc_pa.append(1.)
					if line[7] == '0.0':
						nn+=1
						if int(line[0]) > 0:
							acc.append(0.)
							if line[9] == 'R_e':
								acc_Re.append(0.)
							if line[9] == 'V_e':
								acc_Ve.append(0.)
							if line[9] == 'inc':
								acc_inc.append(0.)
							if line[9] == 'T_p':
								acc_Tp.append(0.)
							if line[9] == 'pa':
								acc_pa.append(0.)
					i=0
					x=0
					lasttime_hrs=''
					lasttime_mins=''
					for i in range(len(line[8])):
						if x == 1:
							lasttime_mins+=line[8][i]
						if line[8][i] == ':':
							x=1
						if x == 0:
							lasttime_hrs+=line[8][i]
					lasttime=float(lasttime_hrs)*3600.+float(lasttime_mins)*60.
				
				else:
					scale=[]
					acc.append(1.)
	#print '{} points'.format(len(this_Re))
	last_pa_diff=0.
	for i in range(len(this_pa)):
		this_pa[i]-=last_pa_diff
		#print i,this_pa[i]
		while this_pa[i] > 90.:
			this_pa[i]-=180.
			last_pa_diff+=180.
		while this_pa[i] < -90.:
			this_pa[i]+=180.
			last_pa_diff-=180.
	
	
	last_pa_diff=0.
	for i in range(len(pa)):
		pa[i]-=last_pa_diff
		while pa[i] > 90.:
			pa[i]-=180.
			last_pa_diff+=180.
		while pa[i] < -90.:
			pa[i]+=180.
			last_pa_diff-=180.
			
	if 'consolidated' not in inp_file:
		cutoff=np.array(nums)[np.where(np.array(chi2s) == min(np.array(chi2s)))][0]-1
	else:
		cutoff=-1
	return nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff

def read_ages(age_file):
	data=ascii.read(age_file)
	anums=np.array(data['Num'])
	age=np.array(data['Age'])*1000.
	mass=np.array(data['Mass'])
	omg_init=np.array(data['omg_init'])
	L_bol=np.array(data['L_bol'])
	L_app=np.array(data['L_app'])
	R_avg=np.array(data['R_avg'])
	R_p=np.array(data['R_p'])
	T_e=np.array(data['T_e'])
	T_avg=np.array(data['T_avg'])
	lg_p=np.array(data['log(g_p)'])
	lg_e=np.array(data['log(g_e)'])
	lg_avg=np.array(data['log(g_avg)'])
	return anums,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg
def plotting(nums,anums,chi2s,R_e,V_e,inc,T_p,pa,acc,model_dir,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg,run_ages,cutoff,out_file):
	nums=np.array(nums)
	chi2s=np.array(chi2s)
	R_e=np.array(R_e)
	V_e=np.array(V_e)
	inc=np.array(inc)
	T_p=np.array(T_p)
	pa=np.array(pa)
	acc=np.array(acc)
	chi2s=chi2s[np.where(nums > cutoff)]
	R_e=R_e[np.where(nums > cutoff)]
	V_e=V_e[np.where(nums > cutoff)]
	inc=inc[np.where(nums > cutoff)]
	T_p=T_p[np.where(nums > cutoff)]
	pa=pa[np.where(nums > cutoff)]
	acc=acc[np.where(nums > cutoff)]
	if run_ages:
		anums=np.array(anums)
		age=np.array(age)
		mass=np.array(mass)
		omg_init=np.array(omg_init)
		L_bol=np.array(L_bol)
		L_app=np.array(L_app)
		R_avg=np.array(R_avg)
		R_p=np.array(R_p)
		T_e=np.array(T_e)
		T_avg=np.array(T_avg)
		lg_p=np.array(lg_p)
		lg_e=np.array(lg_e)
		lg_avg=np.array(lg_avg)
		age=age[np.where(anums > cutoff)]
		mass=mass[np.where(anums > cutoff)]
		omg_init=omg_init[np.where(anums > cutoff)]
		L_bol=L_bol[np.where(anums > cutoff)]
		L_app=L_app[np.where(anums > cutoff)]
		R_avg=R_avg[np.where(anums > cutoff)]
		R_p=R_p[np.where(anums > cutoff)]
		T_e=T_e[np.where(anums > cutoff)]
		T_avg=T_avg[np.where(anums > cutoff)]
		lg_p=lg_p[np.where(anums > cutoff)]
		lg_e=lg_e[np.where(anums > cutoff)]
		lg_avg=lg_avg[np.where(anums > cutoff)]
	nums=nums[np.where(nums > cutoff)]
	
	old_R_e=R_e
	nums=nums[np.where(acc == 1.0)]
	chi2s=chi2s[np.where(acc == 1.0)]
	R_e=R_e[np.where(acc == 1.0)]
	V_e=V_e[np.where(acc == 1.0)]
	inc=inc[np.where(acc == 1.0)]
	T_p=T_p[np.where(acc == 1.0)]
	pa=pa[np.where(acc == 1.0)]
	
	min_chi2=np.amin(chi2s)
	min_chi2_iters=nums[np.where(chi2s == min_chi2)]
	print 'Minimum chi2 is {} and is the chi2 for these iterations: {}'.format(min_chi2,min_chi2_iters)
	open(out_file,'w').write('Minimum chi2 is {} and is the chi2 for these iterations: {}'.format(min_chi2,min_chi2_iters))
	print 'Parameters associated with minimum chi2 are: \n\t{} R_sun, {} km/s, {} deg, {} K, and {} deg'.format(R_e[np.where(chi2s == min_chi2)][0],V_e[np.where(chi2s == min_chi2)][0],inc[np.where(chi2s == min_chi2)][0],T_p[np.where(chi2s == min_chi2)][0],pa[np.where(chi2s == min_chi2)][0])
	open(out_file,'a').write('\nParameters associated with minimum chi2 are: \n\t{} R_sun, {} km/s, {} deg, {} K, and {} deg'.format(R_e[np.where(chi2s == min_chi2)][0],V_e[np.where(chi2s == min_chi2)][0],inc[np.where(chi2s == min_chi2)][0],T_p[np.where(chi2s == min_chi2)][0],pa[np.where(chi2s == min_chi2)][0]))
	
	n_acc=len(R_e)
	
	print '{} points accepted out of the {} runs made after the cutoff, {}'.format(n_acc,len(old_R_e),cutoff+1)
	open(out_file,'a').write('\n{} points accepted out of the {} runs made after the cutoff, {}'.format(n_acc,len(old_R_e),cutoff+1))
	
	for i in range(len(pa)):
		if pa[i] < np.median(pa)-90.:
			pa[i]+=180.
	
	binnum=20
	'''=======================
		Equatorial Radius
	======================='''
	plt.plot(R_e)
	plt.ylabel('Equatorial Radius (R_sun)')
	plt.xlabel('n')
	plt.savefig(model_dir+'plots/conv_re.pdf')
	plt.close()
	
	plt.hist(R_e,bins=binnum)
	plt.xlabel('Equatorial Radius (R_sun)')
	plt.savefig(model_dir+'plots/re_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(R_e)
	sort_re=R_e[this_sort]
	re_68=[]
	re_95=[]
	for i in range(len(sort_re)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			re_68.append(sort_re[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			re_95.append(sort_re[i])
	med_re=np.median(sort_re)
	#print 'R_e: -{}, -{}, {}, +{}, +{}'.format(med_re-min(re_95),med_re-min(re_68),med_re,max(re_68)-med_re,max(re_95)-med_re)
	#open(out_file,'a').write('\nR_e: -{}, -{}, {}, +{}, +{}'.format(med_re-min(re_95),med_re-min(re_68),med_re,max(re_68)-med_re,max(re_95)-med_re))
	print 'R_e: {} || {} - {} || {} - {}'.format(med_re,min(re_68),max(re_68),min(re_95),max(re_95))
	open(out_file,'a').write('R_e: {} || {} - {} || {} - {}'.format(med_re,min(re_68),max(re_68),min(re_95),max(re_95)))
	
	'''==========================
		Equatorial Velocity
	=========================='''
	plt.plot(V_e)
	plt.ylabel('Equatorial Rotation Velocity (km/s)')
	plt.xlabel('n')
	plt.savefig(model_dir+'plots/conv_vel.pdf')
	plt.close()
	
	plt.hist(V_e,bins=binnum)
	plt.xlabel('Equatorial Rotation Velocity (km/s)')
	plt.savefig(model_dir+'plots/vel_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(V_e)
	sort_vel=V_e[this_sort]
	vel_68=[]
	vel_95=[]
	for i in range(len(sort_vel)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			vel_68.append(sort_vel[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			vel_95.append(sort_vel[i])
	med_vel=np.median(sort_vel)
	#print 'V_e: -{}, -{}, {}, +{}, +{}'.format(med_vel-min(vel_95),med_vel-min(vel_68),med_vel,max(vel_68)-med_vel,max(vel_95)-med_vel)
	#open(out_file,'a').write('\nV_e: -{}, -{}, {}, +{}, +{}'.format(med_vel-min(vel_95),med_vel-min(vel_68),med_vel,max(vel_68)-med_vel,max(vel_95)-med_vel))
	print 'V_e: {} || {} - {} || {} - {}'.format(med_vel,min(vel_68),max(vel_68),min(vel_95),max(vel_95))
	open(out_file,'a').write('\nV_e: {} || {} - {} || {} - {}'.format(med_vel,min(vel_68),max(vel_68),min(vel_95),max(vel_95)))
	
	'''=====================
		Inclination
	====================='''
	plt.plot(inc)
	plt.ylabel('Inclination (deg)')
	plt.xlabel('n')
	plt.savefig(model_dir+'plots/conv_inc.pdf')
	plt.close()
	
	plt.hist(inc,bins=binnum)
	plt.xlabel('Inclination (deg)')
	plt.savefig(model_dir+'plots/inc_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(inc)
	sort_inc=inc[this_sort]
	inc_68=[]
	inc_95=[]
	for i in range(len(sort_inc)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			inc_68.append(sort_inc[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			inc_95.append(sort_inc[i])
	med_inc=np.median(sort_inc)
	#print 'inc: -{}, -{}, {}, +{}, +{}'.format(med_inc-min(inc_95),med_inc-min(inc_68),med_inc,max(inc_68)-med_inc,max(inc_95)-med_inc)
	#open(out_file,'a').write('\ninc: -{}, -{}, {}, +{}, +{}'.format(med_inc-min(inc_95),med_inc-min(inc_68),med_inc,max(inc_68)-med_inc,max(inc_95)-med_inc))
	print 'inc: {} || {} - {} || {} - {}'.format(med_inc,min(inc_68),max(inc_68),min(inc_95),max(inc_95))
	open(out_file,'a').write('\ninc: {} || {} - {} || {} - {}'.format(med_inc,min(inc_68),max(inc_68),min(inc_95),max(inc_95)))
	
	'''=====================
		vsini
	====================='''
	vsini=V_e*sin(inc*np.pi/180.)
	
	plt.hist(vsini,bins=binnum)
	plt.xlabel('vsini (km/s)')
	plt.savefig(model_dir+'plots/vsini_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(vsini)
	sort_vsini=vsini[this_sort]
	vsini_68=[]
	vsini_95=[]
	for i in range(len(sort_vsini)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			vsini_68.append(sort_vsini[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			vsini_95.append(sort_vsini[i])
	med_vsini=np.median(sort_vsini)
	#print 'vsini: -{}, -{}, {}, +{}, +{}'.format(med_vsini-min(vsini_95),med_vsini-min(vsini_68),med_vsini,max(vsini_68)-med_vsini,max(vsini_95)-med_vsini)
	#open(out_file,'a').write('\nvsini: -{}, -{}, {}, +{}, +{}'.format(med_vsini-min(vsini_95),med_vsini-min(vsini_68),med_vsini,max(vsini_68)-med_vsini,max(vsini_95)-med_vsini))
	print 'vsini: {} || {} - {} || {} - {}'.format(med_vsini,min(vsini_68),max(vsini_68),min(vsini_95),max(vsini_95))
	open(out_file,'a').write('\nvsini: {} || {} - {} || {} - {}'.format(med_vsini,min(vsini_68),max(vsini_68),min(vsini_95),max(vsini_95)))
		
	'''======================
		Polar Temperature
	======================'''
	plt.plot(T_p)
	plt.ylabel('Polar Temperature (K)')
	plt.xlabel('n')
	plt.savefig(model_dir+'plots/conv_tp.pdf')
	plt.close()
	
	plt.hist(T_p,bins=binnum)
	plt.xlabel('Polar Temperature (K)')
	plt.savefig(model_dir+'plots/tp_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(T_p)
	sort_tp=T_p[this_sort]
	tp_68=[]
	tp_95=[]
	for i in range(len(sort_tp)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			tp_68.append(sort_tp[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			tp_95.append(sort_tp[i])
	med_tp=np.median(sort_tp)
	#print 'T_p: -{}, -{}, {}, +{}, +{}'.format(med_tp-min(tp_95),med_tp-min(tp_68),med_tp,max(tp_68)-med_tp,max(tp_95)-med_tp)
	#open(out_file,'a').write('\nT_p: -{}, -{}, {}, +{}, +{}'.format(med_tp-min(tp_95),med_tp-min(tp_68),med_tp,max(tp_68)-med_tp,max(tp_95)-med_tp))
	print 'T_p: {} || {} - {} || {} - {}'.format(med_tp,min(tp_68),max(tp_68),min(tp_95),max(tp_95))
	open(out_file,'a').write('\nT_p: {} || {} - {} || {} - {}'.format(med_tp,min(tp_68),max(tp_68),min(tp_95),max(tp_95)))
	
	'''===================
		Position Angle
	==================='''
	plt.plot(pa)
	plt.ylabel('Position Angle (deg)')
	plt.xlabel('n')
	plt.savefig(model_dir+'plots/conv_pa.pdf')
	plt.close()
	
	plt.hist(pa,bins=binnum)
	plt.xlabel('Position Angle (deg)')
	plt.savefig(model_dir+'plots/pa_hist.pdf')
	plt.close()
	
	this_sort=np.argsort(pa)
	sort_pa=pa[this_sort]
	pa_68=[]
	pa_95=[]
	for i in range(len(sort_pa)):
		if i+1 > n_acc/6. and i+1 < n_acc*5./6.:
			pa_68.append(sort_pa[i])
		if i+1 > n_acc/40. and i+1 < n_acc*39./40.:
			pa_95.append(sort_pa[i])
	med_pa=np.median(sort_pa)
	#print 'pa: -{}, -{}, {}, +{}, +{}'.format(med_pa-min(pa_95),med_pa-min(pa_68),med_pa,max(pa_68)-med_pa,max(pa_95)-med_pa)
	#open(out_file,'a').write('\npa: -{}, -{}, {}, +{}, +{}'.format(med_pa-min(pa_95),med_pa-min(pa_68),med_pa,max(pa_68)-med_pa,max(pa_95)-med_pa))
	print 'pa: {} || {} - {} || {} - {}'.format(med_pa,min(pa_68),max(pa_68),min(pa_95),max(pa_95))
	open(out_file,'a').write('\npa: {} || {} - {} || {} - {}'.format(med_pa,min(pa_68),max(pa_68),min(pa_95),max(pa_95)))
	if run_ages:
		'''========================
					Extras
		========================'''
		if len(age)+1 != n_acc:
			print '-------------------------'
			print '*****'
			print "Warning. Extras aren't finished yet. {} out of {} done ({}%)".format(len(age),n_acc,round(100.*float(len(age))/float(n_acc),1))
			print '*****'
			open(out_file,'a').write('\n-------------------------')
			open(out_file,'a').write('\n*****')
			open(out_file,'a').write("\nWarning. Extras aren't finished yet. {} out of {} done ({}%)".format(len(age),n_acc,round(100.*float(len(age))/float(n_acc),1)))
			open(out_file,'a').write('\n*****')
		print '-------------------------'
		open(out_file,'a').write('\n-------------------------')
		'''-------
			Age
		-------'''
		num_errors=nums[np.where(age == 0.0)]
		print "The age couldn't be computed in {} instances at the following Nums: {}".format(len(num_errors),num_errors)
		print '-------------------------'
		open(out_file,'a').write("\nThe age couldn't be computed in {} instances at the following Nums: {}".format(len(num_errors),num_errors))
		open(out_file,'a').write('\n-------------------------')
		age=age[np.where(age != 0.0)]
		plt.hist(age,bins=binnum)
		plt.xlabel('Age (Myr)')
		plt.savefig(model_dir+'plots/age_hist.pdf')
		plt.close()
		this_sort=np.argsort(age)
		sort_age=age[this_sort]
		age_68=[]
		age_95=[]
		for i in range(len(sort_age)):
			if i+1 > len(age)/6. and i+1 < len(age)*5./6.:
				age_68.append(sort_age[i])
			if i+1 > len(age)/40. and i+1 < len(age)*39./40.:
				age_95.append(sort_age[i])
		med_age=np.median(sort_age)
		#print 'age: -{}, -{}, {}, +{}, +{}'.format(med_age-min(age_95),med_age-min(age_68),med_age,max(age_68)-med_age,max(age_95)-med_age)
		#open(out_file,'a').write('\nage: -{}, -{}, {}, +{}, +{}'.format(med_age-min(age_95),med_age-min(age_68),med_age,max(age_68)-med_age,max(age_95)-med_age))
		print 'age: {} || {} - {} || {} - {}'.format(med_age,min(age_68),max(age_68),min(age_95),max(age_95))
		open(out_file,'a').write('\nage: {} || {} - {} || {} - {}'.format(med_age,min(age_68),max(age_68),min(age_95),max(age_95)))
		'''-------
			Mass
		-------'''
		mass=mass[np.where(mass != 0.0)]
		plt.hist(mass,bins=binnum)
		plt.xlabel('Mass (M_sun)')
		plt.savefig(model_dir+'plots/mass_hist.pdf')
		plt.close()
		this_sort=np.argsort(mass)
		sort_mass=mass[this_sort]
		mass_68=[]
		mass_95=[]
		for i in range(len(sort_mass)):
			if i+1 > len(mass)/6. and i+1 < len(mass)*5./6.:
				mass_68.append(sort_mass[i])
			if i+1 > len(mass)/40. and i+1 < len(mass)*39./40.:
				mass_95.append(sort_mass[i])
		med_mass=np.median(sort_mass)
		#print 'mass: -{}, -{}, {}, +{}, +{}'.format(med_mass-min(mass_95),med_mass-min(mass_68),med_mass,max(mass_68)-med_mass,max(mass_95)-med_mass)
		#open(out_file,'a').write('\nmass: -{}, -{}, {}, +{}, +{}'.format(med_mass-min(mass_95),med_mass-min(mass_68),med_mass,max(mass_68)-med_mass,max(mass_95)-med_mass))
		print 'mass: {} || {} - {} || {} - {}'.format(med_mass,min(mass_68),max(mass_68),min(mass_95),max(mass_95))
		open(out_file,'a').write('\nmass: {} || {} - {} || {} - {}'.format(med_mass,min(mass_68),max(mass_68),min(mass_95),max(mass_95)))
		'''------------------
			Initial Omega
		------------------'''
		omg_init=omg_init[np.where(omg_init != 0.0)]
		plt.hist(omg_init,bins=binnum)
		plt.xlabel('Initial Omega')
		plt.savefig(model_dir+'plots/omg_init_hist.pdf')
		plt.close()
		this_sort=np.argsort(omg_init)
		sort_omg_init=omg_init[this_sort]
		omg_init_68=[]
		omg_init_95=[]
		for i in range(len(sort_omg_init)):
			if i+1 > len(omg_init)/6. and i+1 < len(omg_init)*5./6.:
				omg_init_68.append(sort_omg_init[i])
			if i+1 > len(omg_init)/40. and i+1 < len(omg_init)*39./40.:
				omg_init_95.append(sort_omg_init[i])
		med_omg_init=np.median(sort_omg_init)
		#print 'omg_init: -{}, -{}, {}, +{}, +{}'.format(med_omg_init-min(omg_init_95),med_omg_init-min(omg_init_68),med_omg_init,max(omg_init_68)-med_omg_init,max(omg_init_95)-med_omg_init)
		#open(out_file,'a').write('\nomg_init: -{}, -{}, {}, +{}, +{}'.format(med_omg_init-min(omg_init_95),med_omg_init-min(omg_init_68),med_omg_init,max(omg_init_68)-med_omg_init,max(omg_init_95)-med_omg_init))
		print 'omg_init: {} || {} - {} || {} - {}'.format(med_omg_init,min(omg_init_68),max(omg_init_68),min(omg_init_95),max(omg_init_95))
		open(out_file,'a').write('\nomg_init: {} || {} - {} || {} - {}'.format(med_omg_init,min(omg_init_68),max(omg_init_68),min(omg_init_95),max(omg_init_95)))
		'''------------------
			Bolometric Luminosity
		------------------'''
		plt.hist(L_bol,bins=binnum)
		plt.xlabel('Bolometric Luminosity (L_sun)')
		plt.savefig(model_dir+'plots/L_bol_hist.pdf')
		plt.close()
		this_sort=np.argsort(L_bol)
		sort_L_bol=L_bol[this_sort]
		L_bol_68=[]
		L_bol_95=[]
		for i in range(len(sort_L_bol)):
			if i+1 > len(L_bol)/6. and i+1 < len(L_bol)*5./6.:
				L_bol_68.append(sort_L_bol[i])
			if i+1 > len(L_bol)/40. and i+1 < len(L_bol)*39./40.:
				L_bol_95.append(sort_L_bol[i])
		med_L_bol=np.median(sort_L_bol)
		#print 'L_bol: -{}, -{}, {}, +{}, +{}'.format(med_L_bol-min(L_bol_95),med_L_bol-min(L_bol_68),med_L_bol,max(L_bol_68)-med_L_bol,max(L_bol_95)-med_L_bol)
		#open(out_file,'a').write('\nL_bol: -{}, -{}, {}, +{}, +{}'.format(med_L_bol-min(L_bol_95),med_L_bol-min(L_bol_68),med_L_bol,max(L_bol_68)-med_L_bol,max(L_bol_95)-med_L_bol))
		print 'L_bol: {} || {} - {} || {} - {}'.format(med_L_bol,min(L_bol_68),max(L_bol_68),min(L_bol_95),max(L_bol_95))
		open(out_file,'a').write('\nL_bol: {} || {} - {} || {} - {}'.format(med_L_bol,min(L_bol_68),max(L_bol_68),min(L_bol_95),max(L_bol_95)))
		'''------------------
			Apparent Luminosity
		------------------'''
		plt.hist(L_app,bins=binnum)
		plt.xlabel('Apparent Luminosity (L_sun)')
		plt.savefig(model_dir+'plots/L_app_hist.pdf')
		plt.close()
		this_sort=np.argsort(L_app)
		sort_L_app=L_app[this_sort]
		L_app_68=[]
		L_app_95=[]
		for i in range(len(sort_L_app)):
			if i+1 > len(L_app)/6. and i+1 < len(L_app)*5./6.:
				L_app_68.append(sort_L_app[i])
			if i+1 > len(L_app)/40. and i+1 < len(L_app)*39./40.:
				L_app_95.append(sort_L_app[i])
		med_L_app=np.median(sort_L_app)
		#print 'L_app: -{}, -{}, {}, +{}, +{}'.format(med_L_app-min(L_app_95),med_L_app-min(L_app_68),med_L_app,max(L_app_68)-med_L_app,max(L_app_95)-med_L_app)
		#open(out_file,'a').write('\nL_app: -{}, -{}, {}, +{}, +{}'.format(med_L_app-min(L_app_95),med_L_app-min(L_app_68),med_L_app,max(L_app_68)-med_L_app,max(L_app_95)-med_L_app))
		print 'L_app: {} || {} - {} || {} - {}'.format(med_L_app,min(L_app_68),max(L_app_68),min(L_app_95),max(L_app_95))
		open(out_file,'a').write('\nL_app: {} || {} - {} || {} - {}'.format(med_L_app,min(L_app_68),max(L_app_68),min(L_app_95),max(L_app_95)))
		'''------------------
			Average Radius
		------------------'''
		plt.hist(R_avg,bins=binnum)
		plt.xlabel('Average Radius (R_sun)')
		plt.savefig(model_dir+'plots/R_avg_hist.pdf')
		plt.close()
		this_sort=np.argsort(R_avg)
		sort_R_avg=R_avg[this_sort]
		R_avg_68=[]
		R_avg_95=[]
		for i in range(len(sort_R_avg)):
			if i+1 > len(R_avg)/6. and i+1 < len(R_avg)*5./6.:
				R_avg_68.append(sort_R_avg[i])
			if i+1 > len(R_avg)/40. and i+1 < len(R_avg)*39./40.:
				R_avg_95.append(sort_R_avg[i])
		med_R_avg=np.median(sort_R_avg)
		#print 'R_avg: -{}, -{}, {}, +{}, +{}'.format(med_R_avg-min(R_avg_95),med_R_avg-min(R_avg_68),med_R_avg,max(R_avg_68)-med_R_avg,max(R_avg_95)-med_R_avg)
		#open(out_file,'a').write('\nR_avg: -{}, -{}, {}, +{}, +{}'.format(med_R_avg-min(R_avg_95),med_R_avg-min(R_avg_68),med_R_avg,max(R_avg_68)-med_R_avg,max(R_avg_95)-med_R_avg))
		print 'R_avg: {} || {} - {} || {} - {}'.format(med_R_avg,min(R_avg_68),max(R_avg_68),min(R_avg_95),max(R_avg_95))
		open(out_file,'a').write('\nR_avg: {} || {} - {} || {} - {}'.format(med_R_avg,min(R_avg_68),max(R_avg_68),min(R_avg_95),max(R_avg_95)))
		'''------------------
			Polar Radius
		------------------'''
		plt.hist(R_p,bins=binnum)
		plt.xlabel('Polar Radius (R_sun)')
		plt.savefig(model_dir+'plots/R_p_hist.pdf')
		plt.close()
		this_sort=np.argsort(R_p)
		sort_R_p=R_p[this_sort]
		R_p_68=[]
		R_p_95=[]
		for i in range(len(sort_R_p)):
			if i+1 > len(R_p)/6. and i+1 < len(R_p)*5./6.:
				R_p_68.append(sort_R_p[i])
			if i+1 > len(R_p)/40. and i+1 < len(R_p)*39./40.:
				R_p_95.append(sort_R_p[i])
		med_R_p=np.median(sort_R_p)
		#print 'R_p: -{}, -{}, {}, +{}, +{}'.format(med_R_p-min(R_p_95),med_R_p-min(R_p_68),med_R_p,max(R_p_68)-med_R_p,max(R_p_95)-med_R_p)
		#open(out_file,'a').write('\nR_p: -{}, -{}, {}, +{}, +{}'.format(med_R_p-min(R_p_95),med_R_p-min(R_p_68),med_R_p,max(R_p_68)-med_R_p,max(R_p_95)-med_R_p))
		print 'R_p: {} || {} - {} || {} - {}'.format(med_R_p,min(R_p_68),max(R_p_68),min(R_p_95),max(R_p_95))
		open(out_file,'a').write('\nR_p: {} || {} - {} || {} - {}'.format(med_R_p,min(R_p_68),max(R_p_68),min(R_p_95),max(R_p_95)))
		'''------------------
			Equatorial Temperature
		------------------'''
		plt.hist(T_e,bins=binnum)
		plt.xlabel('Equatorial Temperature (K)')
		plt.savefig(model_dir+'plots/T_e_hist.pdf')
		plt.close()
		this_sort=np.argsort(T_e)
		sort_T_e=T_e[this_sort]
		T_e_68=[]
		T_e_95=[]
		for i in range(len(sort_T_e)):
			if i+1 > len(T_e)/6. and i+1 < len(T_e)*5./6.:
				T_e_68.append(sort_T_e[i])
			if i+1 > len(T_e)/40. and i+1 < len(T_e)*39./40.:
				T_e_95.append(sort_T_e[i])
		med_T_e=np.median(sort_T_e)
		#print 'T_e: -{}, -{}, {}, +{}, +{}'.format(med_T_e-min(T_e_95),med_T_e-min(T_e_68),med_T_e,max(T_e_68)-med_T_e,max(T_e_95)-med_T_e)
		#open(out_file,'a').write('\nT_e: -{}, -{}, {}, +{}, +{}'.format(med_T_e-min(T_e_95),med_T_e-min(T_e_68),med_T_e,max(T_e_68)-med_T_e,max(T_e_95)-med_T_e))
		print 'T_e: {} || {} - {} || {} - {}'.format(med_T_e,min(T_e_68),max(T_e_68),min(T_e_95),max(T_e_95))
		open(out_file,'a').write('\nT_e: {} || {} - {} || {} - {}'.format(med_T_e,min(T_e_68),max(T_e_68),min(T_e_95),max(T_e_95)))
		'''------------------
			Average Temperature
		------------------'''
		plt.hist(T_avg,bins=binnum)
		plt.xlabel('Average Temperature (K)')
		plt.savefig(model_dir+'plots/T_avg_hist.pdf')
		plt.close()
		this_sort=np.argsort(T_avg)
		sort_T_avg=T_avg[this_sort]
		T_avg_68=[]
		T_avg_95=[]
		for i in range(len(sort_T_avg)):
			if i+1 > len(T_avg)/6. and i+1 < len(T_avg)*5./6.:
				T_avg_68.append(sort_T_avg[i])
			if i+1 > len(T_avg)/40. and i+1 < len(T_avg)*39./40.:
				T_avg_95.append(sort_T_avg[i])
		med_T_avg=np.median(sort_T_avg)
		#print 'T_avg: -{}, -{}, {}, +{}, +{}'.format(med_T_avg-min(T_avg_95),med_T_avg-min(T_avg_68),med_T_avg,max(T_avg_68)-med_T_avg,max(T_avg_95)-med_T_avg)
		#open(out_file,'a').write('\nT_avg: -{}, -{}, {}, +{}, +{}'.format(med_T_avg-min(T_avg_95),med_T_avg-min(T_avg_68),med_T_avg,max(T_avg_68)-med_T_avg,max(T_avg_95)-med_T_avg))
		print 'T_avg: {} || {} - {} || {} - {}'.format(med_T_avg,min(T_avg_68),max(T_avg_68),min(T_avg_95),max(T_avg_95))
		open(out_file,'a').write('\nT_avg: {} || {} - {} || {} - {}'.format(med_T_avg,min(T_avg_68),max(T_avg_68),min(T_avg_95),max(T_avg_95)))
		'''------------------
			log(Polar Surface Gravity)
		------------------'''
		plt.hist(lg_p,bins=binnum)
		plt.xlabel('log(Polar Surface Gravity) (log(cm/s^2))')
		plt.savefig(model_dir+'plots/lg_p_hist.pdf')
		plt.close()
		this_sort=np.argsort(lg_p)
		sort_lg_p=lg_p[this_sort]
		lg_p_68=[]
		lg_p_95=[]
		for i in range(len(sort_lg_p)):
			if i+1 > len(lg_p)/6. and i+1 < len(lg_p)*5./6.:
				lg_p_68.append(sort_lg_p[i])
			if i+1 > len(lg_p)/40. and i+1 < len(lg_p)*39./40.:
				lg_p_95.append(sort_lg_p[i])
		med_lg_p=np.median(sort_lg_p)
		#print 'lg_p: -{}, -{}, {}, +{}, +{}'.format(med_lg_p-min(lg_p_95),med_lg_p-min(lg_p_68),med_lg_p,max(lg_p_68)-med_lg_p,max(lg_p_95)-med_lg_p)
		#open(out_file,'a').write('\nlg_p: -{}, -{}, {}, +{}, +{}'.format(med_lg_p-min(lg_p_95),med_lg_p-min(lg_p_68),med_lg_p,max(lg_p_68)-med_lg_p,max(lg_p_95)-med_lg_p))
		print 'lg_p: {} || {} - {} || {} - {}'.format(med_lg_p,min(lg_p_68),max(lg_p_68),min(lg_p_95),max(lg_p_95))
		open(out_file,'a').write('\nlg_p: {} || {} - {} || {} - {}'.format(med_lg_p,min(lg_p_68),max(lg_p_68),min(lg_p_95),max(lg_p_95)))
		'''------------------
			log(Equatorial Surface Gravity)
		------------------'''
		plt.hist(lg_e,bins=binnum)
		plt.xlabel('log(Equatorial Surface Gravity) (log(cm/s^2))')
		plt.savefig(model_dir+'plots/lg_e_hist.pdf')
		plt.close()
		this_sort=np.argsort(lg_e)
		sort_lg_e=lg_e[this_sort]
		lg_e_68=[]
		lg_e_95=[]
		for i in range(len(sort_lg_e)):
			if i+1 > len(lg_e)/6. and i+1 < len(lg_e)*5./6.:
				lg_e_68.append(sort_lg_e[i])
			if i+1 > len(lg_e)/40. and i+1 < len(lg_e)*39./40.:
				lg_e_95.append(sort_lg_e[i])
		med_lg_e=np.median(sort_lg_e)
		#print 'lg_e: -{}, -{}, {}, +{}, +{}'.format(med_lg_e-min(lg_e_95),med_lg_e-min(lg_e_68),med_lg_e,max(lg_e_68)-med_lg_e,max(lg_e_95)-med_lg_e)
		#open(out_file,'a').write('\nlg_e: -{}, -{}, {}, +{}, +{}'.format(med_lg_e-min(lg_e_95),med_lg_e-min(lg_e_68),med_lg_e,max(lg_e_68)-med_lg_e,max(lg_e_95)-med_lg_e))
		print 'lg_e: {} || {} - {} || {} - {}'.format(med_lg_e,min(lg_e_68),max(lg_e_68),min(lg_e_95),max(lg_e_95))
		open(out_file,'a').write('\nlg_e: {} || {} - {} || {} - {}'.format(med_lg_e,min(lg_e_68),max(lg_e_68),min(lg_e_95),max(lg_e_95)))
		'''------------------
			log(Average Surface Gravity)
		------------------'''
		plt.hist(lg_avg,bins=binnum)
		plt.xlabel('log(Average Surface Gravity) (log(cm/s^2))')
		plt.savefig(model_dir+'plots/lg_avg_hist.pdf')
		plt.close()
		this_sort=np.argsort(lg_avg)
		sort_lg_avg=lg_avg[this_sort]
		lg_avg_68=[]
		lg_avg_95=[]
		for i in range(len(sort_lg_avg)):
			if i+1 > len(lg_avg)/6. and i+1 < len(lg_avg)*5./6.:
				lg_avg_68.append(sort_lg_avg[i])
			if i+1 > len(lg_avg)/40. and i+1 < len(lg_avg)*39./40.:
				lg_avg_95.append(sort_lg_avg[i])
		med_lg_avg=np.median(sort_lg_avg)
		#print 'lg_avg: -{}, -{}, {}, +{}, +{}'.format(med_lg_avg-min(lg_avg_95),med_lg_avg-min(lg_avg_68),med_lg_avg,max(lg_avg_68)-med_lg_avg,max(lg_avg_95)-med_lg_avg)
		#open(out_file,'a').write('\nlg_avg: -{}, -{}, {}, +{}, +{}'.format(med_lg_avg-min(lg_avg_95),med_lg_avg-min(lg_avg_68),med_lg_avg,max(lg_avg_68)-med_lg_avg,max(lg_avg_95)-med_lg_avg))
		print 'lg_avg: {} || {} - {} || {} - {}'.format(med_lg_avg,min(lg_avg_68),max(lg_avg_68),min(lg_avg_95),max(lg_avg_95))
		open(out_file,'a').write('\nlg_avg: {} || {} - {} || {} - {}'.format(med_lg_avg,min(lg_avg_68),max(lg_avg_68),min(lg_avg_95),max(lg_avg_95)))
	'''=================================================
			Postage Stamps
	================================================='''
	the_data=dict()
	the_data['R_e']=R_e
	the_data['V_e']=V_e
	the_data['inc']=inc
	the_data['T_p']=T_p
	the_data['pa']=pa
	
	df = pd.DataFrame(the_data,columns=['R_e','V_e','inc','T_p','pa'])
	
	fig = plt.figure()
	ax = fig.add_subplot(2,1,1)
	scatter_matrix(df, ax=ax,alpha=0.2, figsize=(6, 6), diagonal='hist')
	plt.savefig(model_dir+'plots/correlations.jpg')
	plt.close()
	'''=================================================
			Autocorrelations
	================================================='''
	#corr_re=np.correlate(R_e,R_e,mode='full')
	#corr_vel=np.correlate(V_e,V_e,mode='full')
	#corr_inc=np.correlate(inc,inc,mode='full')
	#corr_tp=np.correlate(T_p,T_p,mode='full')
	#corr_pa=np.correlate(pa,pa,mode='full')
	corr_re=autocorr(R_e)
	corr_vel=autocorr(V_e)
	corr_inc=autocorr(inc)
	corr_tp=autocorr(T_p)
	corr_pa=autocorr(pa)
	
	corr_re/=np.amax(corr_re)
	corr_vel/=np.amax(corr_vel)
	corr_inc/=np.amax(corr_inc)
	corr_tp/=np.amax(corr_tp)
	corr_pa/=np.amax(corr_pa)
	
	plt.plot(corr_re,label='R_e')
	plt.plot(corr_vel,label='V_e')
	plt.plot(corr_inc,label='inc')
	plt.plot(corr_tp,label='T_p')
	plt.plot(corr_pa,label='pa')
	plt.legend()
	plt.savefig(model_dir+'plots/autocorr.pdf')
	plt.close()
	
def autocorr(x):
	x-=np.mean(x)
	result = np.correlate(x, x, mode='full')
	return result[result.size/2:]

if __name__=="__main__":
	main()




