''' Equations taken from Lee et al. (2006)  '''
import numpy as np
''' Constants  '''
uo     = 1.2566370614e-6 # Permeability of free space (m kg s^-2 A^-2)
gH     = 267.513e6       # Proton gyromagnetic ratio (rad s^-1 T^-1)
gN     = -27.126e6       # Nitrogen gyromagnetic ratio (rad s^-1 T^-1)
planck = 6.62606957e-34 # Planck's constant (m^2kg/s)
r_hn   = 1.02e-10        # HN internuclear distance (m)
ddN    = 160.0e-6        # ppm
mg_rat = -9.868831583   # Ratio between H and N 

def tau_c(h_field,two_eta_xy):
  ''' Field values ''' 
  n_field = h_field/mg_rat # this is in Hz
  omega_n   = n_field*(2.0*np.pi) # this is in radians per second
  print 'omega_n = %f'%omega_n
  B0      = (h_field*2.0*np.pi)/gH
  print 'B0 = %f'%B0
  ''' Value of 2eta_xy i.e. R-beta - R-alpha'''
  theta     = 17.0
  theta_rad = np.deg2rad(theta)
  print 'theta_rad = %f'%theta_rad
  ''' Equation 4 '''
  rho = (uo*gH*gN*planck)/(16.0*np.square(np.pi)*np.sqrt(2.0)*np.power(r_hn,3))
  print 'rho = %f'%rho
  ''' Equation 5 '''
  delta_n = (gN*B0*ddN)/(3.0*np.sqrt(2.0))
  print 'delta = %f'%delta_n 
  cos_term = (3.0*(np.square(np.cos(theta_rad)))-1)*(2.0*rho*delta_n)
  x = two_eta_xy/(3.0*(np.square(np.cos(theta_rad)))-1)*(2.0*rho*delta_n)
  print 'x = %f'%x
  ''' Constants for cubic  '''
  a   = (cos_term *1.6 *np.square(omega_n))
  b   = (-two_eta_xy*np.square(omega_n))
  c   = (cos_term*2.8)
  d   = (-two_eta_xy)
  '''Solving cubic equation '''
  coefficients=[a,b,c,d]
  tauc = np.roots(coefficients)
  print 'The tumbling time is %.3f ns'%(np.real(tauc)[0] * 1e9)
  return (np.real(tauc)[0] * 1e9)

if __name__ == "__main__":
  h_field = float(raw_input('What is the field in Hz? '))
  two_eta_xy = float(raw_input('What is the value of Rb-Ra in Hz? '))
  tau_c(h_field,two_eta_xy)
