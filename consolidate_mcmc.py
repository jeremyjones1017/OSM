import read_mcmc as rm
import numpy as np

#stars=['HD110411','HD192640']
#stars=['HD31295','HD110411','HD192640']

grav_darks=['vZ','ELR']
stars=['HD31295','HD110411','HD125162','HD192640']
#stars=['HD125162']
#stars=['HD192640']
#stars=['HD110411']
#stars=['HD31295']

def main(star,grav_dark):
	#grav_dark='vZ'
	#grav_dark='ELR'
	
	star_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'
	if grav_dark == 'vZ':
		models=['mcmc_thread_1','mcmc_thread_2','mcmc_thread_3']
		out_file=star_dir+'consolidated_vZ/'+star+'.mcmc'
	if grav_dark == 'ELR':
		models=['mcmc_thread_4','mcmc_thread_5','mcmc_thread_6']
		out_file=star_dir+'consolidated_ELR/'+star+'.mcmc'
	open(out_file,'w').write('\nNum\tChi^2\tR_e\tV_e\tinc\tT_p\tpa\tacc\tTime\tParam_Changed\trate_100_Re\trate_100_Ve\trate_100_inc\trate_100_Tp\trate_100_pa\trate_Re\trate_Ve\trate_inc\trate_Tp\trate_pa\tscale_Re\tscale_Ve\tscale_inc\tscale_Tp\tscale_pa')

	run_ages=False

	k=0
	for i in models:
		model=i
		model_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+model+'/'
		inp_file=model_dir+star+'.mcmc'
		
		print '--------------------------------'
		print star, model
		print '--------------------------------'
		
		nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff=rm.read(inp_file)
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
		nums=nums[np.where(nums > cutoff)]
		
		old_R_e=R_e
		nums=nums[np.where(acc == 1.0)]
		chi2s=chi2s[np.where(acc == 1.0)]
		R_e=R_e[np.where(acc == 1.0)]
		V_e=V_e[np.where(acc == 1.0)]
		inc=inc[np.where(acc == 1.0)]
		T_p=T_p[np.where(acc == 1.0)]
		pa=pa[np.where(acc == 1.0)]
		acc=acc[np.where(acc == 1.0)]
	
		
		min_chi2=np.amin(chi2s)
		min_chi2_iters=nums[np.where(chi2s == min_chi2)]
		print 'Minimum chi2 is {} and is the chi2 for these iterations: {}'.format(min_chi2,min_chi2_iters)
		print 'Parameters associated with minimum chi2 are: \n\t{} R_sun, {} km/s, {} deg, {} K, and {} deg'.format(R_e[np.where(chi2s == min_chi2)][0],V_e[np.where(chi2s == min_chi2)][0],inc[np.where(chi2s == min_chi2)][0],T_p[np.where(chi2s == min_chi2)][0],pa[np.where(chi2s == min_chi2)][0])
		
		n_acc=len(R_e)
		
		print '{} points accepted out of the {} runs made after the cutoff, {}'.format(n_acc,len(old_R_e),cutoff+1)
		
		
		for i in range(len(pa)):
			if pa[i] < np.median(pa)-90.:
				pa[i]+=180.
				
		for j in range(len(nums)):
			open(out_file,'a').write('\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(j+k,chi2s[j],R_e[j],V_e[j],inc[j],T_p[j],pa[j],acc[j],'--1','--2','--3','--4','--5','--6','--7','--8','--9','--10','--11','--12','--13','--14','--15','--16','--17'))
		k+=len(nums)

	if grav_dark == 'vZ':
		model='consolidated_vZ'
	if grav_dark == 'ELR':
		model='consolidated_ELR'
	
	model_dir='C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/'+star+'/'+model+'/'
	inp_file=model_dir+star+'.mcmc'
	age_file=model_dir+star+'.ages'
	out_file=model_dir+star+'.results'
	print '--------------------------------'
	print star, model
	print '--------------------------------'
	nums,chi2s,R_e,V_e,inc,T_p,pa,acc,acc_Re,acc_Ve,acc_inc,acc_Tp,acc_pa,scale,nn,cutoff=rm.read(inp_file)
	if run_ages:
		anums,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg=rm.read_ages(age_file)
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
	rm.plotting(nums,anums,chi2s,R_e,V_e,inc,T_p,pa,acc,model_dir,age,mass,omg_init,L_bol,L_app,R_avg,R_p,T_e,T_avg,lg_p,lg_e,lg_avg,run_ages,cutoff,out_file)
	
	
if __name__=="__main__":
	#main()
	for j in grav_darks:
		for i in stars:
			main(i,j)




