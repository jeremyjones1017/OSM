using Interpolations
using PyPlot

#Constants
NG=6.67384e-8 #Newton's Gravity in cm^3/g/s^2
R_sun=6.955e10 #Solar Radius in cm
M_sun=1.988435e33 #Solar Mass in g
L_sun=3.839e33 #Solar Luminosity in erg/s
pc=3.08567758e18 #1 parsec in cm
sigma_SB=5.6704e-5 #Stefan-Boltzmann constant in erg/cm^2/s/K^4
h=6.626e-27 #Planck's constant in cm^2*g/s
c=3e10 #Speed of light, cm/s
k=1.381e-16 #Boltzmann constant erg/K

#=
	LIST OF FUNCTIONS TO TRANSCRIBE		[The functions to transcribe in here that they call]
	osm									[calc_beta,sort_hull_results,read_phoenix,calc_vis,calc_phot,calc_Lbol,age_mass]
	read_phoenix						[do_phx_integrate]
	extract_phoenix_full				[read_this_phoenix,read_this_phoenix_ftp]
	extract_phoenix_phot				[read_this_phoenix,read_this_phoenix_ftp]
	extract_phoenix_vis					[read_this_phoenix,read_this_phoenix_ftp]
	sort_hull_results					[N/A]
	extract								[extract_phoenix_vis]
	calc_beta							[ftht,froche]
	ftht								[N/A]
	froche								[N/A]
	age_mass							[match,amoeba]
	match								[N/A]
	read_this_phoenix					[do_phx_integrate]
	read_this_phoenix_ftp				[do_phx_integrate]
	amoeba								[N/A]
	calc_vis							[extract,plot_ellipse,plot_vis]
	calc_phot							[extract_phoenix_phot,plot_phot]
	calc_Lbol							[extract_phoenix_full]
	read_cwlzpf							[N/A]
	get_phoenix_wave					[N/A]
	get_complext_trf					[N/A]
	do_phx_integrate					[N/A]
	read_input							[N/A]
	plot_ellipse						[N/A]
	plot_vis							[N/A]
	plot_phot							[extract_phoenix_full]
=#

function read_filters(use_filts,filt_dir,cwl,wav)
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
	filt_dict=Dict()
	for i in 1:length(use_filts)
		this_wav_tran=readdlm(filt_dir*use_filts[i]*".txt")
		this_wav=this_wav_tran[:,1]/1e8
		this_tran=this_wav_tran[:,2]
		int_tran=zeros(wav)
		
		itp=interpolate((this_wav,),this_tran,Gridded(Linear()))	#This interpolates the saved filter curve to the appropriate wavelength scale
		for j in 1:length(wav)
			if wav[j] >= this_wav[1] && wav[j] <= this_wav[end]
				int_tran[j]=itp[wav[j]]
			end
		end
		abs_wav_minus_cwl=abs(wav - cwl[use_filts[i]])
		cwl_tran=int_tran[abs_wav_minus_cwl .== minimum(abs_wav_minus_cwl)]
		int_tran=int_tran/cwl_tran[1]
		filt_dict[use_filts[i]]=int_tran
	end
	return filt_dict
end
function fwhm(wave,transmission)
	#=Determines the full width half max of the supplied transmission curve
	Inputs:
	wave
		An array of wavelengths
	transmission
		The fractional transmission of the filter curve at wavelengths 'wave'
	Output:
	fwhm
		The full width half max of the transmission curve
	=#
	#The following test to see if the transmission curve is a top hat function rather 
	#than a normal transmission curve
	test_trans=transmission[transmission .> 0.01]	#An array of all the nonzero points in transmission
	test_trans=test_trans[test_trans .< 0.99]		#An array of all the nonzero/non-one points in the transmission
	if length(test_trans) == 0	#If the tranmsission array is only zeros and ones
		ones_wave=wave[transmission .> 0.01]
		fwhm_high=maximum(ones_wave)
		fwhm_low=minimum(ones_wave)
		fwhm=fwhm_high-fwhm_low	#The fwhm is the difference between where the ones end and where they begin
	else	#If the transmission array isn't a top hat function
		close_to_one=minimum((transmission-1).^2)	#I'm not sure why I did it this way... I should double check this later
		max_wave=wave[(transmission-1).^2 .== close_to_one]
		if isodd(length(max_wave))
			max_wave=max_wave[convert(Int64,length(max_wave)/2+0.5)]
		else
			max_wave=max_wave[convert(Int64,length(max_wave)/2)]
		end
		wave_high=wave[wave .> max_wave]
		trans_high=transmission[wave .> max_wave]
		close_to_half=minimum((trans_high-0.5).^2)
		fwhm_high=wave_high[(trans_high-0.5).^2 .== close_to_half]
		wave_low=wave[wave .< max_wave]
		trans_low=transmission[wave .< max_wave]
		close_to_half=minimum((trans_low-0.5).^2)
		fwhm_low=wave_low[(trans_low-0.5).^2 .== close_to_half]
		fwhm=fwhm_high-fwhm_low
	end
	return fwhm[1]
end
function unitrange(res)
	#=Outputs an array with values ranging from 0 to 1 with a number of elements given by the input.
	Input:
	res
		The number of elements you want in the array
	Output:
	return
		An array that ranges from 0 to 1 with res elements
	=#
	return collect(linspace(0,1,res))
end
function read_cwlzpf(cwlzpf_file)
	data=readdlm(cwlzpf_file)
	filt=data[:,1]
	shift!(filt)		#Removes 'title' from array
	cwl=data[:,2]
	shift!(cwl)
	zpf=data[:,3]
	shift!(zpf)
	cwl_dict=Dict()
	zpf_dict=Dict()
	for i in 1:length(filt)
		cwl_dict[filt[i]]=float(cwl[i])
		zpf_dict[filt[i]]=float(zpf[i])
	end
	return cwl_dict,zpf_dict
end
function cart2sphere(xxx,yyy,zzz,inc,pa)
	"""Converts the input cartesian coordinates into spherical coordinates (adjusting for inclination and position angle of the star)
	Inputs:
	xxx
		The x coordinate to be converted
	yyy
		The y coordinate to be converted
	zzz
		The z coordinate to be converted
	inc
		The inclination of the model star
	pa
		The position angle of the model star
	
	Outputs:
	r
		The associated r coordinate
	tht
		The associated tht coordinate
	phi
		The associated phi coordinate
	"""
	#Rotating by the position angle
	xx=xxx*cos(pa)+yyy*sin(pa)
	yy=-xxx*sin(pa)+yyy*cos(pa)
	zz=zzz
	#Rotating by the inclination
	x=xx
	y=yy*sin(inc)+zz*cos(inc)
	z=-yy*cos(inc)+zz*sin(inc)
	#Converting to spherical coordinates
	r=sqrt(x^2+y^2+z^2)
	phi=acos(z/sqrt(z^2+x^2))
	tht=acos(y/sqrt(z^2+x^2+y^2))
	return r,tht,phi
end
function store_mesa(fil)
	"""Collects the mass track from a given file and gets the relevant info from it and stores it as a list ready
		to be packaged into mesa_dict
	Inputs:
	fil
		The input file name
		
	Outputs:	The following are in the returned list
	age
		An array with the ages (in Gyr) of the mass track
	teff
		An array with the log(T_eff/K) of the mass track
	lum
		An array with the log(L/L_sun) of the mass track
	rad
		An array with the log(R/R_sun) of the mass track
	r_p
		An array with the log(R_p/R_sun) of the mass track
	vel
		An array with the Equatorial Velocity in km/s of the mass track
	wnow
		An array with the Angular rotation rate/critical at the current age of the mass track
	"""
	data=readdlm(fil)
	age=data[:,2]
	age=Array{Float64}(age[6:end])		#Age in yr
	age=age*1e-9						#Age in Gyr
	teff=data[:,4]
	teff=Array{Float64}(teff[6:end])	#log(T_eff/K)
	lum=data[:,5]
	lum=Array{Float64}(lum[6:end])		#log(L/L_sun)
	rad=data[:,6]
	rad=Array{Float64}(rad[6:end])		#log(R/R_sun)
	vel=data[:,15]
	vel=Array{Float64}(vel[6:end])		#Equatorial Velocity in km/s
	wnow=data[:,16]
	wnow=Array{Float64}(wnow[6:end])	#Angular rotation rate/critical at the current age
	return [age,teff,lum,rad,vel,wnow]
end
function read_mesa(masses,omegas,mesa_dir,mesa_use_Z)
	"""Creates mesa_dict - the dictionary of mass tracks stored for the given masses/omegas combo
	
	Inputs:
	masses
		An array of strings representing the masses available for the given mesa_use_Z
	omegas
		An array of strings representing the omegas available for the given mesa_use_Z
	mesa_dir
		The directory in which the mesa output files are stored
	mesa_use_Z
		The internal metallicity to be used
	
	Outputs:
	mesa_dict
		The dictionary of mass tracks stored for the given masses/omegas combo	
	"""
	mesa_dict=Dict{String,Any}()
	for i in 1:length(masses)
		for j in 1:length(omegas)
			inp_file=mesa_dir*mesa_use_Z*"_M"*masses[i]*"_w"*omegas[j]
			#if i == 1 && j == 1
			#	mesa_dict=Dict{String,Any}(inp_file => store_mesa(inp_file))
			#else
			#	println(mesa_dict)
			mesa_dict[inp_file]=store_mesa(inp_file)
			#end
		end
	end
	return mesa_dict
end
function read_vis(vis_inp)
	data=readdlm(vis_inp)
	wl=Array{Float64}(data[:,1])
	wlerr=Array{Float64}(data[:,2])
	vis=Array{Float64}(data[:,3])
	vis_err=Array{Float64}(data[:,4])
	u_m=Array{Float64}(data[:,5])
	v_m=Array{Float64}(data[:,6])
	u_l=Array{Float64}(data[:,7])
	v_l=Array{Float64}(data[:,8])
	cal=Array{String}(data[:,9])
	#Convert wavelengths from um to m
	wlerr*=1e-6
	wl*=1e-6
	#Sort the data by wavelength
	the_sort=sortperm(wl)
	vis=vis[the_sort]
	vis_err=vis_err[the_sort]
	u_m=u_m[the_sort]
	v_m=v_m[the_sort]
	u_l=u_l[the_sort]
	v_l=v_l[the_sort]
	cal=cal[the_sort]
	wlerr=wlerr[the_sort]
	wl=wl[the_sort]
	return wl,wlerr,vis,vis_err,u_m,v_m,u_l,v_l,cal
end
function read_phot(phot_inp)
	phot_data=Dict{String,Array{Float64}}()
	data=readdlm(phot_inp)
	use_filts=Array{String}(data[:,1])
	phot=Array{Float64}(data[:,2])
	phot_err=Array{Float64}(data[:,3])
	for i in 1:length(use_filts)
		#println(use_filts[i],' ',phot[i],' ',phot_err[i])
		phot_data[use_filts[i]]=[phot[i],phot_err[i]]
	end
	return phot_data,use_filts
end
function read_input(filename)
	data_dict=Dict{String,Any}()
	data=readdlm(filename)
	println(data)
	#=
	with open(filename,'r') as input:
		input_reader=csv.reader(input,delimiter='|',skipinitialspace=True)
		input_reader.next()
		for line in input_reader:
			for i in range(len(line)):
				line[i] = re.sub('\t', '', line[i])
			data_dict[line[0]]=line[1]
			
	star=data_dict['Star']
	model=data_dict['Model']
	star_dir=data_dict['Star Directory']+star+'/'
	model_dir=star_dir+model+'/'
	confirm_file=model_dir+star+'.confirm'
	open(confirm_file,'a').write('Time is {}'.format(time.ctime(time.time())))
	open(confirm_file,'w').write('\n Input file {} has been read and is running. \n The following flags have been set:'.format(filename))
	if data_dict['Calc Vis'] == 'Y': open(confirm_file,'a').write('\n Visibilities will be calculated')
	if data_dict['Calc Phot'] == 'Y': open(confirm_file,'a').write('\n Photometry will be calculated')
	if data_dict['Calc Lum'] == 'Y': open(confirm_file,'a').write('\n Luminosity will be calculated')
	if data_dict['Calc Age'] == 'Y': open(confirm_file,'a').write('\n Age will be calculated')
	if data_dict['GPU Accel'] == 'Y': open(confirm_file,'a').write('\n GPU acceleration will be used')
	if data_dict['Verbose'] == 'Y': open(confirm_file,'a').write('\n Verbose mode will be used')
	if data_dict['Gravity Darkening'] == 'vZ': open(confirm_file,'a').write('\n The vZ gravity darkening law will be used')
	if data_dict['Gravity Darkening'] == 'ELR': open(confirm_file,'a').write('\n The ELR gravity darkening law will be used')

	return data_dict
	=#
end

read_input("C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/HD192640/HD192640_win.input")



#=
phot_inp="C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/HD192640/HD192640.phot"
phot_data,use_filts=read_phot(phot_inp)

for i in phot_data
	println(i[1],' ',phot_data[i[1]])
end
vis_inp="C:/Users/Jeremy/Dropbox/Programing/Astars/Stars/HD192640/HD192640.vis"
read_vis(vis_inp)

mesa_dir="C:/Users/Jeremy/Dropbox/Programing/Astars/MESA/History_Files/"
mesa_use_Z="Z0.0111"	#[M/H]=-0.14
fmasses=Vector(linspace(1.0,3.3,24))
fomegas=Vector(linspace(0.0,0.9,10))
masses=Vector{String}(length(fmasses))
omegas=Vector{String}(length(fomegas))
for i in 1:length(fmasses)
	masses[i]=string(fmasses[i])
end
for i in 1:length(fomegas)
	omegas[i]=string(fomegas[i])
end

read_mesa(masses,omegas,mesa_dir,mesa_use_Z)

cwl,zpf=read_cwlzpf("C:/Users/Jeremy/Dropbox/Programing/Astars/Band_Passes/cwlzpf.txt")
use_filts=[]
for i in cwl
	push!(use_filts,i[1])
end
#use_filts=["2massH","2massJ","2massK","cousinsI","cousinsR","F1565","F1965","F2365","F2740","johnsonB","johnsonH","johnsonI","johnsonJ","johnsonK","johnsonR","johnsonR","johnsonU","johnsonV","sloanr","sloanu","stromgrenb","stromgrenu","stromgrenv","stromgreny","wes15N","wes15W","wes18","wes22","wes25","wes33"]
#use_filts=["2massJ","stromgrenb","stromgreny","wes15N","wes15W","wes18","wes22","wes25","wes33"]
phx_wav=Vector(linspace(500e-8,25999e-8,25500))
filt_dict=read_filters(use_filts,"C:/Users/Jeremy/Dropbox/Programing/Astars/Band_Passes/",cwl,phx_wav)
for i in filt_dict
	println(i[1]," ",sum(filt_dict[i[1]]))
end
=#