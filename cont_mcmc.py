import OSMlib as osm
import read_mcmc as rm
import mcmc_osm as mc
import os
import numpy as np
import time
import random as rn

def main():
	total_time_start=time.time()
	bashCommand = "cls"
	os.system(bashCommand)

	input_file='/nfs/morgan/users/jones/Dropbox/Python/Astars/Stars/HD125162/mcmc_thread_4/HD125162_mcmc_thread_4.input'
	input_dict=osm.read_input(input_file)
	star=input_dict['Star']
	model=input_dict['Model']
	star_dir=input_dict['Star Directory']+star+'/'
	model_dir=star_dir+model+'/'
	inp_file=model_dir+star+'.mcmc'
	phx_dir=input_dict['Atmo Directory']
	filt_dir=input_dict['Filter Directory']
	
	nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff=rm.read(inp_file)
	
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
		
	end_model_at=10000
	#end_model_at=6985
	for i in range(len(nums)):
		if nums[i] < cutoff and acc[i] == 1.:
			end_model_at+=1
	
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
	
	base_chi2=chi2s[-1]
	m=float(input_dict['Mass'])
	#m=2.06219230493	Used for HD 106951 mcmc_thread_1
	#m=2.098			Used for HD 106951 mcmc_thread_2/3
	#m=1.8				Used for HD 110411 (all)
	
	beta=0.
	dist=1000./float(input_dict['Parallax'])
	vsini=float(input_dict['vsini'])
	vsini_err=float(input_dict['vsini_err'])
	
	
	lock_Re=False
	lock_vel=False
	lock_inc=False
	lock_Tp=False
	lock_pa=False
	
	the_params=[]
	free_params=[]
	if lock_Re == False:
		free_params.append(0)
		the_params.append('R_e')
	if lock_vel == False:
		free_params.append(1)
		the_params.append('V_e')
	if lock_inc == False:
		free_params.append(2)
		the_params.append('inc')
	if lock_Tp == False:
		free_params.append(3)
		the_params.append('T_p')
	if lock_pa == False:
		free_params.append(4)
		the_params.append('pa')
	
	r=[]
	for i in range(len(nums)):
		r.append([R_e[i],V_e[i],inc[i]*np.pi/180.,T_p[i],pa[i]*np.pi/180.+np.pi/2.])
	g_scale=1.
	
	empty_phx_dict=dict()
	first_chi2,phx_dict,g_points,extras=osm.osm(r[0],[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,empty_phx_dict,colat_len,phi_len,mode])
	if 'v' in mode:
		print 'Increasing the scale so that the number of points is ~1000. This will be the base chi2'
		g_scale=np.sqrt(1000./g_points)
		first_chi2,phx_dict,g_points,extras=osm.osm(r[0],[first_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode])
	
	n=nums[-1]+1
	mc.mcmc(r,n,nn,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,total_time_start,end_model_at,free_params,base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode,vsini,vsini_err,inp_file,the_params)
	

if __name__=="__main__":
	main()
