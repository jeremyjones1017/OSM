import OSMlib as osm
import os
import numpy as np
import time

def main():
	total_time_start=time.time()
	bashCommand = "cls"
	os.system(bashCommand)
	
	input_dict=osm.read_input('C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/HD192640/HD192640_win.input')
	
	star=input_dict['Star']
	model=input_dict['Model']
	star_dir=input_dict['Star Directory']+star+'/'
	model_dir=star_dir+model+'/'
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
	#if input_dict['Doppler'] == 'Y': mode+='d'
	
	wl,wlerr,vis,vis_err,u_m,v_m,u_l,v_l,cal=osm.read_vis(vis_inp)
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
	
	R_e=float(input_dict['Equatorial Radius'])
	V_e=float(input_dict['Equatorial Velocity'])
	Inc=float(input_dict['Inclination'])
	T_p=float(input_dict['Polar Temperature'])
	PA=float(input_dict['Position Angle'])
	vsini=float(input_dict['vsini'])
	vsini_err=float(input_dict['vsini_err'])
	
	m=float(input_dict['Mass'])
	beta=0.		#Holdover from old version that I haven't gotten rid of. This beta is meaningless

	dist=1000./float(input_dict['Parallax'])

	g_scale=1.
	
	r1=[R_e,V_e,Inc*np.pi/180.,T_p,PA*np.pi/180.+np.pi/2.]
	r2=[R_e,V_e,Inc*np.pi/180.,T_p+250.,PA*np.pi/180.+np.pi/2.]
	r3=[R_e,V_e,Inc*np.pi/180.,T_p+500.,PA*np.pi/180.+np.pi/2.]
	r4=[R_e,V_e,Inc*np.pi/180.,T_p+750.,PA*np.pi/180.+np.pi/2.]
	r5=[R_e,V_e,Inc*np.pi/180.,T_p+1000.,PA*np.pi/180.+np.pi/2.]
	
	phx_dict=dict()
	chi2,phx_dict,g_points,extras=osm.osm(r1,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	g_scale=np.sqrt(1000./g_points)

	#if input_dict['Do Plot'] == 'Y': mode+='P'

	chi2,phx_dict,g_points,extras=osm.osm(r1,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	chi2,phx_dict,g_points,extras=osm.osm(r2,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	chi2,phx_dict,g_points,extras=osm.osm(r3,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	chi2,phx_dict,g_points,extras=osm.osm(r4,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	chi2,phx_dict,g_points,extras=osm.osm(r5,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,u_m,v_m,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,cwl,phx_dict,colat_len,phi_len,cal,star,model,model_dir,mode])
	
	
	total_time_finish=time.time()
	total_time=total_time_finish-total_time_start
	ttm=int(total_time/60.)
	tts=total_time-ttm*60.
	print 'Done: {} m {} s elapsed'.format(ttm,tts)
	
if __name__=="__main__":
	main()