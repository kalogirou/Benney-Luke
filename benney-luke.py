
# Solve Benney-Luke equations

from firedrake import *
from ufl import *
from parameters import *
from mesh import *
from IC_etaR import *
from solvers import *
from energy import *
from exact import *

# A bug in Firedrake kernel optimiser currently means that this code fails with optimisations enabled, so we turn it off
op2.init()
parameters["coffee"]["O2"] = False

# Parameter and experimental values
(is_exact, is_linear, order_time, solvers_print) = get_cases();
(H0s, hh0, hh1, LLx, LLy, Lc, LL, U, g) = experimental_values();
(T, dt, mu, epsilon, Nx, Ny, Lx, Ly, d, Ts) = parameter_values(H0s, LLx, LLy, Lc, hh1, U, g, is_exact, is_linear);
if is_exact == 'True': (Ampl, kx, ky, k2, omega, speed) = omega_speed(Lx, Ly, mu);

# Create mesh 
(mesh, coords) = create_mesh(Nx, Ny, Lx, Ly);
if is_exact == 'False': coords = stretch_mesh(mesh, coords, Nx, Ny, Lx, d);

# Define functions
V = FunctionSpace(mesh, "CG", 2)

eta0 = Function(V)
phi0 = Function(V)
eta1 = Function(V)
phi1 = Function(V)
eta0_5 = Function(V)
phi0_5 = Function(V)
q0 = Function(V)
q1 = Function(V)
q0_5 = Function(V)
etaR = Function(V)
etaR0_5 = Function(V)
ex_eta = Function(V)
ex_phi = Function(V)

eta = TrialFunction(V)
phi = TrialFunction(V)
q = TrialFunction(V)
gamma = TestFunction(V)

# Initial conditions
if is_exact == 'True':
	(eta0, phi0) = exact_ICs(eta0, phi0, coords.dat.data, kx, ky, Lx, Ly, Ampl, mu, epsilon, speed, is_linear);	# Exact solutions.
elif is_exact == 'False':
	(etaR_expr, h0, h1) = def_etaR(coords.dat.data, H0s, mu, epsilon, hh0, hh1, LL);
	etaR.interpolate(etaR_expr)			# Define etaR
	eta0 = eta0_eq_etaR(eta,eta0,etaR,gamma);	# Solve LinearVariationalProblem eta=etaR for initial eta0

if is_linear == 'True': epsilon = 0.0;			# Set epsilon=0 for linear case

# Weak formulation: order_time = 1 for Euler 1st-order, = 2 for Stormer-Verlet 2nd-order
if order_time == 1:
	phi_solver = solver_phi(phi1, phi1, phi0, eta0, etaR, phi, gamma, dt, mu, epsilon, is_linear, solvers_print);
	(q_solver, q_solver0_5) = solvers_q(q1, q0_5, phi1, phi0_5, q, gamma);
	eta_solver = solver_eta(eta1, eta0, eta0, phi1, phi1, q1, q1, eta, gamma, dt, mu, epsilon, is_linear, solvers_print);		# Give eta0 instead of eta0_5 to cancel with coef. 0.5

elif order_time == 2:
	phi_solver0_5 = solver_phi(phi0_5, phi0_5, phi0, eta0, etaR0_5, phi, gamma, 0.5*dt, mu, epsilon, is_linear, solvers_print);	# Give phi0_5 instead of phi1
	(q_solver1, q_solver0_5) = solvers_q(q1, q0_5, phi1, phi0_5, q, gamma);
	eta_solver1 = solver_eta(eta1, eta1, eta0, phi1, phi0_5, q1, q0_5, eta, gamma, dt, mu, epsilon, is_linear, solvers_print);	# Give eta1 instead of eta0_5
	phi_solver1 = solver_phi(phi1, phi0_5, phi0_5, eta1, etaR0_5, phi, gamma, 0.5*dt, mu, epsilon, is_linear, solvers_print);	# Give eta1 instead of eta0 and phi0_5 instead of phi0

# Write data to files
phi_file = File("phi.pvd")
eta_file = File("eta.pvd")
if is_exact == 'False': etaR_file = File("etaR.pvd")

phi_file << phi0
eta_file << eta0
if is_exact == 'False': etaR_file << etaR

# Time iteration
t = 0
set_E0 = 0

# Initial energy
E0 = assemble( ( 0.5*(1+epsilon*eta0)*abs(grad(phi0))**2 + 0.5*eta0**2 + mu*(inner(grad(q0),grad(phi0))-0.75*q0**2) )*dx )
E_file = open("energy.txt", "w")
E = find_energy(phi0, eta0, q0, mu, epsilon, t, E0, E_file);

# Exact solutions
if is_exact == 'True': 
	(expr_phi, expr_eta) = exact_solutions(coords.dat.data, t, kx, ky, k2, Lx, Ly, omega, Ampl, mu, epsilon, speed, is_linear);
	
	phi_exact = File("phi_ex.pvd")					# Write exact solutions to files
	eta_exact = File("eta_ex.pvd")
	phi_exact << phi0
	eta_exact << eta0

while(t < T-dt):
    	t += dt
    	if is_exact == 'True':
    		print t, assemble((eta0-ex_eta)**2*dx)

		# Update exact solutions in time
		(ex_phi, ex_eta) = update_exact(ex_phi, ex_eta, expr_phi, expr_eta, t, Ly, epsilon, speed, is_linear);

		phi_exact << ex_phi					# Write exact solutions to files
		eta_exact << ex_eta
	elif is_exact == 'False':
		print t, E

		# Update etaR in time
		if t<Ts:
	    		etaR_expr.H0 = h1+(h0-h1)*(Ts-t)/Ts		
    		else:
	    		etaR_expr.H0 = h1

    		etaR.interpolate(etaR_expr)

    		# Evaluate etaR at t_{n+1/2}=t_n+dt/2
    		if t<Ts:
	    		etaR_expr.H0 = h1+(h0-h1)*(Ts-(t+0.5*dt))/Ts

    		etaR0_5.interpolate(etaR_expr)

	# Solve the variational problem
	if order_time == 1:
		(phi0, eta0, q0) = solvers_SE(phi_solver, q_solver, eta_solver, phi0, eta0, q0, phi1, eta1, q1);
	elif order_time == 2:
		(phi0, eta0, q0) = solvers_SV(phi_solver0_5, q_solver0_5, eta_solver1, phi_solver1, q_solver1, phi0, eta0, q0, phi1, eta1, q1);

    	# Monitor energy in time
	if is_exact == 'False':
    		if t>Ts:
			# As soon the sluice gate is released (i.e. t>Ts), calculate E0 again, close E_file and open it again so that only the energy values from this point on are saved
			set_E0 = set_E0 + 1
			if set_E0 == 1:
				E0 = assemble( ( 0.5*(1+epsilon*eta0)*abs(grad(phi0))**2 + 0.5*eta0**2 + mu*(inner(grad(q0),grad(phi0))-0.75*q0**2) )*dx )
				E_file.close()
				E_file = open("energy.txt", "w")

    	E = find_energy(phi0, eta0, q0, mu, epsilon, t, E0, E_file);

	# Write data to files
    	phi_file << phi0
    	eta_file << eta0
	if is_exact == 'False': etaR_file << etaR

E_file.close()
