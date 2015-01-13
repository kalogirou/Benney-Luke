
from math import *


##################################################

def get_cases():

	is_exact = 'False'	# 'True': compare to exact solutions. 'False': sluice gate problem.
	is_linear = 'False'	# 'True': linear Benney-Luke. 'False': nonlinear.
	order_time = 2		# Order of the time integrator: = 1 Symplectic Euler, = 2 Stormer-Verlet.
	solvers_print = {} 	# OR = {'snes_monitor': True,'ksp_monitor': True,'snes_linesearch_monitor': True}

	return (is_exact, is_linear, order_time, solvers_print);

##################################################

def experimental_values():

	H0s = 0.45	# m	 # average water level at rest
	hh0 = 0.43	# m	 # lower water level at rest
	hh1 = 0.9	# m	 # upper water lever at rest
	LLx = 2.0	# m	 # wavetank width
	LLy = 43.63	# m	 # wavetank length
	Lc  = 2.7	# m	 # contraction length
	LL  = 2.63	# m	 # location of sluice gate
	U   = 2.5	# m/s	 # sluice gate release speed
	g   = 9.81	# m/s^2	 # gravity

	return (H0s, hh0, hh1, LLx, LLy, Lc, LL, U, g);

##################################################

def parameter_values(H0s, LLx, LLy, Lc, hh1, U, g, is_exact, is_linear):

	dt = 0.0028
	mu = 0.04
	epsilon = 0.25

	if is_exact == 'True':
		d = 0.0
		Ts = 0.0
		if is_linear == 'True':			# Linear exact problem.
			T = 2.0
			Nx = 200
			Ny = 36
			Lx = 5
			Ly = 1.8
		elif is_linear == 'False':		# Soliton exact problem.
			T = 10.0
			Nx = 1
			Ny = 50
			Lx = 1
			Ly = 10
			mu = 0.01
			epsilon = 0.01 
	elif is_exact == 'False':			# Sluice gate problem.
		T = 25.0
		Nx = 20
		Ny = 800
		Lx = LLx/(H0s/sqrt(mu))		# width (scaled)
		Ly = LLy/(H0s/sqrt(mu))		# length (scaled)
		d = Lc/(H0s/sqrt(mu))		# contraction length (scaled)
		Ts = (hh1/U)/sqrt(H0s/g/mu)	# sluice gate release time

	return (T, dt, mu, epsilon, Nx, Ny, Lx, Ly, d, Ts);

##################################################

def omega_speed(Lx, Ly, mu):

	Ampl = 0.1
	kx = 2
	ky = 4
	k2 = (kx*pi/Lx)**2 + (ky*pi/Ly)**2
	omega = sqrt(k2*(1+2*mu/3.0*k2))/(1+0.5*mu*k2)	# Dispersion relation

	# Necessary only in the nonlinear soliton case
	speed = 1.5 					# corresponding to c=1 in the solution in paper -> c is scaled differently

	return (Ampl, kx, ky, k2, omega, speed);

##################################################

