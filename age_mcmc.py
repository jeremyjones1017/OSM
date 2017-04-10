import OSMlib as osm
import read_mcmc as rm
import mcmc_osm as mc
import os
import os.path
import numpy as np
import time
import random as rn
import csv

def main():
	#start_nums=input('What Num to start with? ')
	total_time_start=time.time()
	bashCommand = "cls"
	os.system(bashCommand)

	input_dict=osm.read_input('C:/Users/Jeremy/Dropbox/Python/Astars/Stars/HD125162/HD125162_win.input')

	star=input_dict['Star']
	model=input_dict['Model']
	star_dir=input_dict['Star Directory']+star+'/'
	model_dir=star_dir+model+'/'
	phx_dir=input_dict['Atmo Directory']
	filt_dir=input_dict['Filter Directory']
	inp_file=model_dir+star+'.mcmc'
	out_file=model_dir+star+'_M-0.14.ages'
	use_Z='Z-0.0'
	
	nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff=rm.read(inp_file,'ages')
	
	vis_inp=star_dir+star+'.vis'	#Visibility input file
	phot_inp=star_dir+star+'.phot'	#Photometry input file
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
	
	base_chi2=chi2s[-1]
	m=float(input_dict['Mass'])
	beta=0.		#Holdover from old version that I haven't gotten rid of. This beta is meaningless
	dist=1000./float(input_dict['Parallax'])
	vsini=float(input_dict['vsini'])
	vsini_err=float(input_dict['vsini_err'])
	
	if os.path.exists(out_file):
		print 'file exists'
		anums=[]
		with open(out_file,'r') as input:
			input_reader=csv.reader(input,delimiter='\t')
			input_reader.next()
			for line in input_reader:
				anums.append(int(line[0]))
		start_nums=max(anums)+1
	else:
		start_nums=0
	
	print 'Num\tAge\tMass\tomg_init\tL_bol\tL_app\tR_avg\tR_p\tT_e\tT_avg\tlog(g_p)\tlog(g_e)\tlog(g_avg)'
	if start_nums == 0:
		open(out_file,'w').write('Num\tAge\tMass\tomg_init\tL_bol\tL_app\tR_avg\tR_p\tT_e\tT_avg\tlog(g_p)\tlog(g_e)\tlog(g_avg)')
	
	phx_dict=dict()
	g_scale=1.
	
	
	nums=np.array(nums)
	chi2s=np.array(chi2s)
	R_e=np.array(R_e)
	V_e=np.array(V_e)
	inc=np.array(inc)
	T_p=np.array(T_p)
	pa=np.array(pa)
	acc=np.array(acc)
	nums=nums[np.where(acc == 1.0)]
	chi2s=chi2s[np.where(acc == 1.0)]
	R_e=R_e[np.where(acc == 1.0)]
	V_e=V_e[np.where(acc == 1.0)]
	inc=inc[np.where(acc == 1.0)]
	T_p=T_p[np.where(acc == 1.0)]
	pa=pa[np.where(acc == 1.0)]
	acc=acc[np.where(acc == 1.0)]
	
	
	i=0
	while nums[i] < start_nums:
		i+=1
	while i < len(nums)-1:
		r=[R_e[i],V_e[i],inc[i]*np.pi/180.,T_p[i],pa[i]*np.pi/180.+np.pi/2.]
		if i % 20 == 0 and i != 0:
			print '------------------------------------------------------------------------------------'
			print '{} ages out of {} calculated ({}%)'.format(i,len(nums),round(100.*float(i+1)/float(len(nums)),1))
			print '------------------------------------------------------------------------------------'
			print 'Num\tAge\tMass\tomg_init\tL_bol\tL_app\tR_avg\tR_p\tT_e\tT_avg\tlog(g_p)\tlog(g_e)\tlog(g_avg)'
			print '------------------------------------------------------------------------------------'
		if i == 0:
			#print nums[i],chi2s[i]
			chi2,phx_dict,g_points,extras=osm.osm(r,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode])
			print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8])
			open(out_file,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8]))
		else:
			if pa[i] != pa[i-1]:
				#print 'pa is what changed, so no effect on these parameters'
				print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8])
				open(out_file,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8]))
			#print nums[i],chi2s[i]
			else:
				chi2,phx_dict,g_points,extras=osm.osm(r,[base_chi2,m,beta,dist,vis,vis_err,phot_data,wl,u_l,v_l,uni_wl,uni_dwl,g_scale,phx_dir,use_Z,use_filts,filt_dict,zpf,phx_dict,colat_len,phi_len,mode])
				print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8])
				open(out_file,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(nums[i],extras[9],extras[10],extras[11],extras[0],extras[1],extras[2],extras[3],extras[4],extras[5],extras[6],extras[7],extras[8]))
		i+=1
		
	

if __name__=="__main__":
	main()