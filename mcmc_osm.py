import OSMlib as osm
import os
import numpy as np
import time
import random as rn

def main():
	total_time_start=time.time()
	bashCommand = "cls"
	os.system(bashCommand)
	
	input_file='/nfs/morgan/users/jones/Dropbox/Python/Astars/Stars/HD31295/mcmc_thread_4/HD31295_mcmc_thread_4.input'
	input_dict=osm.read_input(input_file)
	
	star=input_dict['Star']
	model=input_dict['Model']
	star_dir=input_dict['Star Directory']+star+'/'
	model_dir=star_dir+model+'/'
	rot_out=model_dir+star+'.mcmc'
	phx_dir=input_dict['Atmo Directory']
	filt_dir=input_dict['Filter Directory']
	use_Z='Z-0.0'
	vis_inp=star_dir+star+'.vis'	#Visibility input file
	phot_inp=star_dir+star+'.phot'	#Photometry input file
	#mode='vgo'	#v-calculate visibilities, p-calculate photometry, L - calc Luminosities, r-record results (text and plots), g-use gpu, o-print outputs, a-calculate ages/masses
	mode=''
	if input_dict['Calc Vis'] == 'Y': mode+='v'
	if input_dict['Calc Phot'] == 'Y': mode+='p'
	if input_dict['Calc Lum'] == 'Y': mode+='L'
	if input_dict['Calc Age'] == 'Y': mode+='a'
	if input_dict['GPU Accel'] == 'Y': mode+='g'
	if input_dict['Verbose'] == 'Y': mode+='o'
	if input_dict['Gravity Darkening'] == 'vZ': mode+='z'
	if input_dict['Gravity Darkening'] == 'ELR': mode+='r'
	
	wl,wlerr,vis,vis_err,u_m,v_m,u_l,v_l=osm.read_vis(vis_inp)
	phot_data,use_filts=osm.read_phot(phot_inp)
	cwl,zpf=osm.read_cwlzpf(filt_dir+'cwlzpf.txt')
	if input_dict['Atmo Model'] == 'Phoenix': wav=osm.get_phoenix_wave(phx_dir)
	filt_dict=osm.read_filters(use_filts,filt_dir,cwl,wav)
	
	colat_len=20
	phi_len=30
	
	uni_wl=[]	#What are all the unique wavlengths in this observation
	uni_dwl=[] #The fwhm of the unique wavelengths observed
	for i in range(len(wl)):
		if wl[i] not in uni_wl:
			uni_wl.append(wl[i])
			uni_dwl.append(wlerr[i])
	
	base_chi2=1e8
	m=float(input_dict['Mass'])
	#m=2.06219230493	Used for HD 106951 mcmc_thread_1
	#m=2.098			Used for HD 106951 mcmc_thread_2/3
	#m=1.8				Used for HD 110411 (all)
	
	beta=0.
	dist=1000./float(input_dict['Parallax'])
	
	
	g_scale=1.
	
	R_e=float(input_dict['Equatorial Radius'])
	V_e=float(input_dict['Equatorial Velocity'])
	Inc=float(input_dict['Inclination'])
	T_p=float(input_dict['Polar Temperature'])
	PA=float(input_dict['Position Angle'])
	vsini=float(input_dict['vsini'])
	vsini_err=float(input_dict['vsini_err'])

	end_model_at=10000
	
	r=[]
	r.append([R_e,V_e,Inc*np.pi/180.,T_p,PA*np.pi/180.+np.pi/2.])
	
	g_scale=1.
	
	empty_phx_dict=dict()
	base_chi2,phx_dict,g_points,extras=osm.osm(r[0],[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,empty_phx_dict,colat_len,phi_len,mode])
	if 'v' in mode:
		print 'Increasing the scale so that the number of points is ~1000. This will be the base chi2'
		g_scale=np.sqrt(1000./g_points)
		base_chi2,phx_dict,g_points,extras=osm.osm(r[0],[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode])

	scale=[0.08,15.,3.*np.pi/180.,130.,43.*np.pi/180.] #The initial range for mcmc to search over
	scale=np.array(scale)
	
	lock_Re=False
	lock_vel=False
	lock_inc=False
	lock_Tp=False
	lock_pa=False
	
	the_params=['R_e','V_e','inc','T_p','pa']
	free_params=[]
	if lock_Re == False:
		free_params.append(0)
	if lock_vel == False:
		free_params.append(1)
	if lock_inc == False:
		free_params.append(2)
	if lock_Tp == False:
		free_params.append(3)
	if lock_pa == False:
		free_params.append(4)
	
	n=1
	nn=0
	acc=[1.]
	acc_Re=[1.]
	acc_Ve=[1.]
	acc_inc=[1.]
	acc_Tp=[1.]
	acc_pa=[1.]
	curtime=time.time()
	elapsed=(curtime-total_time_start)/60.
	elapsed_hrs=0
	while elapsed > 60.:
		elapsed-=60.
		elapsed_hrs+=1
	print '\nNum\tChi^2\tR_e\tV_e\tinc\tT_p\tpa\tacc\tacc_all\tTime'
	open(rot_out,'w').write('\nNum\tChi^2\tR_e\tV_e\tinc\tT_p\tpa\tacc\tTime\tParam_Changed\trate_100_Re\trate_100_Ve\trate_100_inc\trate_100_Tp\trate_100_pa\trate_Re\trate_Ve\trate_inc\trate_Tp\trate_pa\tscale_Re\tscale_Ve\tscale_inc\tscale_Tp\tscale_pa')
	print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(0,base_chi2,r[0][0],r[0][1],r[0][2]*180./np.pi,r[0][3],r[0][4]*180./np.pi-90.,'base','--',str(elapsed_hrs)+':'+str(elapsed))
	open(rot_out,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(0,base_chi2,r[0][0],r[0][1],r[0][2]*180./np.pi,r[0][3],r[0][4]*180./np.pi-90.,'base',str(elapsed_hrs)+':'+str(elapsed),'--','--','--','--','--','--','--','--','--','--','--',scale[0],scale[1],scale[2]*180./np.pi,scale[3],scale[4]*180./np.pi))
	mcmc(r,n,nn,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,total_time_start,end_model_at,free_params,base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode,vsini,vsini_err,rot_out,the_params)

def mcmc(r,n,nn,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,total_time_start,end_model_at,free_params,base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode,vsini,vsini_err,rot_out,the_params):
	while sum(acc)+1.<=end_model_at:
		n_acc=int(sum(acc)+1)
		if acc[-1] == 1. and n_acc % 20 == 0:
			print '------------------------------------------------------------------------------------'
			print '{} accepted models out of {} ({}%)'.format(n_acc,end_model_at,round(100.*float(n_acc)/float(end_model_at),2))
			print '------------------------------------------------------------------------------------'
			print 'Num\tChi^2\tR_e\tV_e\tinc\tT_p\tpa\tacc\tacc_all\tTime'
			print '------------------------------------------------------------------------------------'
		if scale[4] > np.pi/2.:
			scale[4] = np.pi/2.
		this_param=rn.choice(free_params)
		this_Re=r[n-nn-1][0]
		this_vel=r[n-nn-1][1]
		this_inc=r[n-nn-1][2]
		this_Tp=r[n-nn-1][3]
		this_pa=r[n-nn-1][4]
		last_vsini = this_vel*np.sin(this_inc)
		if this_param == 0:
			this_Re=rn.gauss(r[n-nn-1][0],scale[0])
		if this_param == 1:
			this_vel=rn.gauss(r[n-nn-1][1],scale[1])
		if this_param == 2:
			this_inc=rn.gauss(r[n-nn-1][2],scale[2])
			while this_inc > np.pi/2. or this_inc < 0.002:
				this_inc=rn.gauss(r[n-nn-1][2],scale[2])
		if this_param == 3:
			this_Tp=rn.gauss(r[n-nn-1][3],scale[3])
		if this_param == 4:
			this_pa=rn.gauss(r[n-nn-1][4],scale[4])
		while this_pa > np.pi:
			this_pa-=np.pi
			#print 'pa decreased to {}'.format(this_pa)
		while this_pa < 0.:
			this_pa+=np.pi
			#print 'pa increased to {}'.format(this_pa)
		
		this_r=[this_Re,this_vel,this_inc,this_Tp,this_pa]
		r.append(this_r)
		this_chi2,phx_dict,g_points,extras=osm.osm(this_r,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode])
		curtime=time.time()
		elapsed=(curtime-total_time_start)/60.
		elapsed_hrs=0
		while elapsed > 60.:
			elapsed-=60.
			elapsed_hrs+=1
		this_vsini = this_vel*np.sin(this_inc)
		#a=np.exp(-0.5*(this_chi2-base_chi2))*gaussian(this_vsini,vsini,vsini_err)
		this_prior=gaussian(this_vsini,vsini,vsini_err)
		last_prior=gaussian(last_vsini,vsini,vsini_err)
		#a=np.exp(-(this_chi2-base_chi2-this_prior+last_prior))		#This ends up being a weak prior on vsini - leads to ridiculous vsini distribution given the observed vsini
		a=np.exp(-(this_chi2-base_chi2))*(this_prior/last_prior)	#This is a stronger prior on vsini.
		if a > 1:
			#print '{} km/s, {} km/s, {} km/s'.format(vsini-2.*vsini_err,this_vsini,vsini+2.*vsini_err)
			base_chi2=this_chi2
			acc.append(1.)
			nn=0
			if this_param == 0:
				acc_Re.append(1.)
			if this_param == 1:
				acc_Ve.append(1.)
			if this_param == 2:
				acc_inc.append(1.)
			if this_param == 3:
				acc_Tp.append(1.)
			if this_param == 4:
				acc_pa.append(1.)
		else:
			if a > rn.random():
				#print '{} km/s, {} km/s, {} km/s'.format(vsini-2.*vsini_err,this_vsini,vsini+2.*vsini_err)
				base_chi2=this_chi2
				acc.append(1.)
				if this_param == 0:
					acc_Re.append(1.)
				if this_param == 1:
					acc_Ve.append(1.)
				if this_param == 2:
					acc_inc.append(1.)
				if this_param == 3:
					acc_Tp.append(1.)
				if this_param == 4:
					acc_pa.append(1.)
				nn=0
			else:
				nn+=1
				acc.append(0.)
				if this_param == 0:
					acc_Re.append(0.)
				if this_param == 1:
					acc_Ve.append(0.)
				if this_param == 2:
					acc_inc.append(0.)
				if this_param == 3:
					acc_Tp.append(0.)
				if this_param == 4:
					acc_pa.append(0.)

		if len(acc_Re) < 100:
			rate_Re_100 = np.average(acc_Re)
		else:
			rate_Re_100=0.
			for i in range(100):
				rate_Re_100+=acc_Re[-1-i]
			rate_Re_100/=100.			
		if len(acc_Ve) < 100:
			rate_Ve_100 = np.average(acc_Ve)
		else:
			rate_Ve_100=0.
			for i in range(100):
				rate_Ve_100+=acc_Ve[-1-i]
			rate_Ve_100/=100.			
		if len(acc_inc) < 100:
			rate_inc_100 = np.average(acc_inc)
		else:
			rate_inc_100=0.
			for i in range(100):
				rate_inc_100+=acc_inc[-1-i]
			rate_inc_100/=100.			
		if len(acc_Tp) < 100:
			rate_Tp_100 = np.average(acc_Tp)
		else:
			rate_Tp_100=0.
			for i in range(100):
				rate_Tp_100+=acc_Tp[-1-i]
			rate_Tp_100/=100.			
		if len(acc_pa) < 100:
			rate_pa_100 = np.average(acc_pa)
		else:
			rate_pa_100=0.
			for i in range(100):
				rate_pa_100+=acc_pa[-1-i]
			rate_pa_100/=100.
		
		if len(acc_Re) > 0:
			rate_Re=np.average(acc_Re)
			if len(acc_Re) % 20 == 0:
				if this_param == 0:
					if rate_Re_100 > 0.5:
						scale[0]*=1.5
					elif rate_Re_100 < 0.1:
						scale[0]*=0.5
		else:
			rate_Re=0.
			rate_Re_100=0.
		if len(acc_Ve) > 0:
			rate_Ve=np.average(acc_Ve)
			if len(acc_Ve) % 20 == 0:
				if this_param == 1:
					if rate_Ve_100 > 0.5:
						scale[1]*=1.5
					elif rate_Ve_100 < 0.1:
						scale[1]*=0.5
		else:
			rate_Ve=0.
			rate_Ve_100=0.
		if len(acc_inc) > 0:
			rate_inc=np.average(acc_inc)
			if len(acc_inc) % 20 == 0:
				if this_param == 2:
					if rate_inc_100 > 0.5:
						scale[2]*=1.5
					elif rate_inc_100 < 0.1:
						scale[2]*=0.5
		else:
			rate_inc=0.
			rate_inc_100=0.
		if len(acc_Tp) > 0:
			rate_Tp=np.average(acc_Tp)
			if len(acc_Tp) % 20 == 0:
				if this_param == 3:
					if rate_Tp_100 > 0.5:
						scale[3]*=1.5
					elif rate_Tp_100 < 0.1:
						scale[3]*=0.5
		else:
			rate_Tp=0.
			rate_Tp_100=0.
		if len(acc_pa) > 0.:
			rate_pa=np.average(acc_pa)
			if len(acc_pa) % 20 == 0:
				if this_param == 4:
					if rate_pa_100 > 0.5:
						scale[4]*=1.5
					elif rate_pa_100 < 0.1:
						scale[4]*=0.5
		else:
			rate_pa=0.
			rate_pa_100=0.
		
		print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(n,base_chi2,r[n-nn][0],r[n-nn][1],r[n-nn][2]*180./np.pi,r[n-nn][3],r[n-nn][4]*180./np.pi-90.,acc[n],np.average(acc),str(elapsed_hrs)+':'+str(elapsed))
		open(rot_out,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(n,base_chi2,r[n-nn][0],r[n-nn][1],r[n-nn][2]*180./np.pi,r[n-nn][3],r[n-nn][4]*180./np.pi-90.,acc[n],str(elapsed_hrs)+':'+str(elapsed),the_params[this_param],rate_Re_100,rate_Ve_100,rate_inc_100,rate_Tp_100,rate_pa_100,rate_Re,rate_Ve,rate_inc,rate_Tp,rate_pa,scale[0],scale[1],scale[2]*180./np.pi,scale[3],scale[4]*180./np.pi))
		n+=1

	
	total_time_finish=time.time()
	total_time=total_time_finish-total_time_start
	ttm=int(total_time/60.)
	tts=total_time-ttm*60.
	print 'Done: {} m {} s elapsed'.format(ttm,tts)

def gaussian(x, mu, sig):
    return np.exp(-(x - mu)**2. / (2 * sig**2.))

if __name__=="__main__":
	main()
