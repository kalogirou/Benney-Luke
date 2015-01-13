
from firedrake import *


##################################################

def find_energy(phi0, eta0, q0, mu, epsilon, t, E0, E_file):

	# Energy calculated at every time step
	E = assemble( ( 0.5*(1+epsilon*eta0)*abs(grad(phi0))**2 + 0.5*eta0**2 + mu*(inner(grad(q0),grad(phi0))-0.75*q0**2) )*dx )
    	time = str(t)
    	Et = str(E)
    	Ed = str(abs(E-E0)/E0)		# Energy difference from initial value E0
    	E_file.write('%-10s %-10s %-10s\n' % (time,Ed,Et))

	return E;

##################################################

