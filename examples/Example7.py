import numpy as np
from enterprise_GWecc import GWecc
import matplotlib.pyplot as plt

year = 365.25*24*3600
day = 24*3600
parsec = 102927125.0
MSun = 4.92703806e-6
ns = 1e-9

def hms_to_rad(hh,mm,ss):
	sgn = np.sign(hh)
	return sgn * (sgn*hh + mm/60 + ss/3600) * np.pi/12

def dms_to_rad(dd,mm,ss):
	sgn = np.sign(dd)
	return sgn * (sgn*dd + mm/60 + ss/3600) * np.pi/180

M = 1e9
q = 1
Pb0 = 1.5	# years
#n0 = 2*np.pi/Pb
Omega = 0
i = 0
t0 = 0
l0 = gamma0 = 0
z = 0

RA_P  = hms_to_rad(  4,37,15.81476)
DEC_P = dms_to_rad(-47,15, 8.6242)
D_P = 156.3	# pc
     
RA_GW = hms_to_rad(  4,0,0)
DEC_GW = dms_to_rad(-45,0,0)
D_GW = 1e9 # pc 

ntoas = 5000
toas = 365.25*np.linspace(0,10,ntoas)	# days

cosmu, Fp, Fx = GWecc.AntennaPattern(RA_GW, DEC_GW, RA_P, DEC_P)
delay = -1000#-D_P*(1-cosmu) / (1+z)

term = GWecc.ResidualsTerms_Earth

for idx,e0 in enumerate([0.1,0.3,0.6]):

	res_px_anl = GWecc.EccentricResiduals_px(M, q,
					                       	 Omega, i,
					                         t0, Pb0, e0, l0, gamma0,
					                       	 D_GW, delay,
					                       	 z,
					                       	 GWecc.ResidualsMethod_Anl,
					                       	 term,
					                         toas)

	res_px_num = GWecc.EccentricResiduals_px(M, q,
					                       	 Omega, i,
					                         t0, Pb0, e0, l0, gamma0,
					                       	 D_GW, delay,
					                       	 z,
					                       	 GWecc.ResidualsMethod_Num,
					                       	 term,
					                         toas)
	
	res_p_anl = np.asarray(res_px_anl[0]) - np.mean(res_px_anl[0])
	res_x_anl = np.asarray(res_px_anl[1]) - np.mean(res_px_anl[1])
	
	res_p_num = np.asarray(res_px_num[0]) - np.mean(res_px_num[0])
	res_x_num = np.asarray(res_px_num[1]) - np.mean(res_px_num[1])
	
	ax = plt.subplot(321+idx*2)
	plt.plot(toas/365.25, res_p_num, label="$s_+$ Analytic", color="b")
	plt.plot(toas/365.25, res_p_anl, "--", label="$s_+$ Numerical", color="g")
	plt.grid()
	plt.xlim([0,10])
	plt.ylabel("$s_{+,\\times}(t_E)$  (ns)", fontsize=14)
	if idx==0:
		plt.legend(loc='upper right', fontsize=12)
	
	if idx==2:
		plt.tick_params(labelsize=12)
		plt.xlabel("$t_E$ (years)", fontsize=14)
	else:
		plt.tick_params(labelbottom=False,labelsize=12)
		
	plt.text(1,0.8*ax.yaxis.get_data_interval()[1],"e="+str(e0),size=15, ha="center", va="center", bbox=dict(boxstyle="round",facecolor='cyan',alpha=0.5))
	
	ax = plt.subplot(321+idx*2+1)
	plt.plot(toas/365.25, res_x_num, label="$s_\\times$ Analytic", color="b")
	plt.plot(toas/365.25, res_x_anl, "--", label="$s_\\times$ Numerical", color="g")
	plt.grid()
	plt.xlim([0,10])
	#plt.ylabel("$s_{+,\\times}(t)$  (ns)", fontsize=14)
	if idx==0:
		plt.legend(loc='upper right', fontsize=12)
	
	if idx==2:
		plt.tick_params(labelsize=12)
		plt.xlabel("$t=t_E$ (years)", fontsize=14)
	else:
		plt.tick_params(labelbottom=False,labelsize=12)
	
	"""
	for jdx,method in enumerate([GWecc.ResidualsMethod_Num, GWecc.ResidualsMethod_Anl]):
		
		

		
		
		
		
		
		if term == GWecc.ResidualsTerms_Earth:
			
			
		
		if idx == 0:
			title = "Earth term" if term==GWecc.ResidualsTerms_Earth else "Pulsar term"
			plt.title(title)
			
		
		
		if idx<2:
			
		else:
			
	"""
plt.show()
