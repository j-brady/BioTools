''' script to calculate Mw from Tau_c values - by Jacob Brady'''
import numpy as np
import pylab as pl

pi   = np.pi		
#eta  = 631e-6 		# viscosity of water at 315K
eta  = 855e-6 		# viscosity of water at 300K
#r_w  = 1.6e-10 		# to 3.2 hydration radius
r_w  = 3.2e-10 		# to 3.2 hydration radius
N_a  = 6.0221413e23	# avagadro's number
#rho  = 1.37e-6 		# m^3g^-1
rho  = 0.73e-6 		# m^3g^-1
Mw   = 16.181e3		# protein molecular weight
k    = 1.3806488e-23 	# boltzmann contant  m^2 kg s^-2 K^-1
#T    = 315		# temp in kelvin
T    = 300		# temp in kelvin
#tau_c= 7.89e-9		# s
#Ubiquitin
tau_c= 3.8e-9		# s

def viscosity(T): # Temperature is in Kelvin
  ''' Viscosity calculation for water '''
  mu_0 =  0.001792        # kg/ms
  T_0  =  273.16          # Kelvin
  a    = -1.94  
  b    = -4.80
  c    =  6.74
  mu = mu_0 * np.exp((a+b*(T_0/T)+c*np.square(T_0/T)))
  return mu

def mw_from_tauc(r_w,eta,T,tau_c,rho):
  equation = (3*tau_c*k*T)/(4*pi*eta)
  r = equation**(1/3.0)
  print('The effective hydrodynamic radius = %.3f A'%(r*1e10))
  m_w = ((r-r_w)**3.0)*4.0*pi*N_a/(3.0*rho)
  print('Molecular weight = %.3f kDa'% (m_w/1000.0))
  return m_w


if __name__ == "__main__":
  # r is the hydrodynamic radius
  #Ubiquitin
  fig = pl.subplot(111)
  fig.set_xlabel(r'$hydration$ $layer$ $(\AA)$',fontsize=20)
  fig.set_ylabel('$molecular$ $weight$ $(kDa)$',fontsize=20)

  r_w_array = np.arange(1.6e-10,3.21e-10,0.01e-10)
  micelle_array = []
  bicelle_array = []
  for i in r_w_array:
    #micelle = mw_from_tauc(r_w=i,eta=viscosity(315),T=315,tau_c=7.89e-9,rho=rho)
    #bicelle = mw_from_tauc(r_w=i,eta=viscosity(315),T=315,tau_c=11.87e-9,rho=rho)
    micelle = mw_from_tauc(r_w=i,eta=viscosity(310),T=310,tau_c=7.17e-9,rho=rho)
    bicelle = mw_from_tauc(r_w=i,eta=viscosity(314),T=314,tau_c=7.17e-9,rho=rho)
    micelle_array.append(micelle)
    bicelle_array.append(bicelle)
  
  micelle_kDa = np.array(micelle_array)/1000.0
  bicelle_kDa = np.array(bicelle_array)/1000.0
  r_w_angstrom= r_w_array*1e10

  #fig.plot(r_w_angstrom,micelle_kDa,c='r',label='$micelle$')
  #fig.plot(r_w_angstrom,bicelle_kDa,c='b',label='$bicelle$')
  fig.plot(r_w_angstrom,micelle_kDa,c='r',label='$37$')
  fig.plot(r_w_angstrom,bicelle_kDa,c='b',label='$41$')
  fig.set_xlim(1.6,3.2)
  leg = pl.legend(prop={'size':20})
  leg.draw_frame(False)
  pl.tight_layout(pad=0.8,h_pad=0.8,w_pad=0.4)
  pl.tick_params(axis='both', which='major', labelsize=16)
  pl.savefig('molecular_weights.pdf')
  ##pl.show()    
  print(viscosity(315))
