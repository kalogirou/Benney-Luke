
from firedrake import *


##################################################

# Variational problem for phi equation
def solver_phi(phi1, phi0_5, phi0, eta0, etaR, phi, gamma, dt, mu, epsilon, is_linear, solvers_print):

	# Note: For the Euler method, give   phi1    instead of	phi0_5.
	# Note: For the SV method, give      0.5*dt	  "	dt,
	#			   and	     etaR0_5	  "	etaR.
	# 	-> If in addition give 	     phi0_5	  "	phi1,	then this is the first half-step of the SV.
	#	-> If in addition give	     phi0_5	  "	phi0,
	#		          and	     eta1  	  "	eta0,	then this is the second half-step of the SV.
	if is_linear == 'True':
		aphi = ( gamma*phi + 0.5*mu*inner(grad(gamma),grad(phi)) )*dx
		Lphi = ( gamma*phi0 + 0.5*mu*inner(grad(gamma),grad(phi0)) - dt*gamma*(eta0-etaR) )*dx

		phi_problem = LinearVariationalProblem(aphi,Lphi,phi1)
		phi_solver = LinearVariationalSolver(phi_problem)
	elif is_linear == 'False':
		Fphi = ( gamma*(phi1-phi0)/dt + 0.5*mu*inner(grad(gamma),grad((phi1-phi0)/dt)) + gamma*(eta0-etaR) + 0.5*epsilon*inner(grad(phi0_5),grad(phi0_5))*gamma )*dx

		phi_problem = NonlinearVariationalProblem(Fphi,phi1)
		phi_solver = NonlinearVariationalSolver(phi_problem,solver_parameters=solvers_print)

	return phi_solver;

# Variational problem for equation q = -2/3*Delta^2(phi)
def solvers_q(q1, q0_5, phi1, phi0_5, q, gamma):

	aq = gamma*q*dx
	Lq = 2.0/3.0*inner(grad(gamma),grad(phi1))*dx
	Lq0_5 = 2/3.0*inner(grad(gamma),grad(phi0_5))*dx

	q_problem = LinearVariationalProblem(aq,Lq,q1)
	q_solver  = LinearVariationalSolver(q_problem)

	q_problem0_5 = LinearVariationalProblem(aq,Lq0_5,q0_5)	
	q_solver0_5  = LinearVariationalSolver(q_problem0_5)

	return (q_solver, q_solver0_5);

# Variational problem for eta equation
def solver_eta(eta1, eta0_5, eta0, phi1, phi0_5, q1, q0_5, eta, gamma, dt, mu, epsilon, is_linear, solvers_print):

	# Note: For the Euler method, give   eta0    instead of	eta0_5,
	#			      and    phi1	  "	phi0_5,
	#			      and    q1		  "	q0_5.
	# Note: For the SV method, give      eta1	  "	eta0_5.
	if is_linear == 'True':
		aeta = ( gamma*eta + 0.5*mu*inner(grad(gamma),grad(eta)) )*dx
		Leta = ( gamma*eta0 + 0.5*mu*inner(grad(gamma),grad(eta0)) + dt*inner(grad(gamma),grad(phi0_5)) + mu*dt*inner(grad(gamma),grad(q0_5)) )*dx

		eta_problem = LinearVariationalProblem(aeta,Leta,eta1)
		eta_solver = LinearVariationalSolver(eta_problem)
	elif is_linear == 'False':
		Feta = ( gamma*(eta1-eta0)/dt + 0.5*mu*inner(grad(gamma),grad((eta1-eta0)/dt)) - 0.5*((1+epsilon*eta0)+(1+epsilon*eta0_5))*inner(grad(gamma),grad(phi0_5)) - mu*inner(grad(gamma),grad(q0_5)) )*dx

		eta_problem = NonlinearVariationalProblem(Feta,eta1)
		eta_solver = NonlinearVariationalSolver(eta_problem,solver_parameters=solvers_print)

	return eta_solver;

##################################################

# 1st-order Symplectic Euler solvers
def solvers_SE(phi_solver, q_solver, eta_solver, phi0, eta0, q0, phi1, eta1, q1):

	phi_solver.solve()
	q_solver.solve()    
	eta_solver.solve()

	eta0.assign(eta1)
	phi0.assign(phi1)
	q0.assign(q1)

	return (phi0, eta0, q0);

##################################################

# 2nd-order Stormer-Verlet solvers
def solvers_SV(phi_solver0_5, q_solver0_5, eta_solver1, phi_solver1, q_solver1, phi0, eta0, q0, phi1, eta1, q1):

    	phi_solver0_5.solve()
    	q_solver0_5.solve()    
    	eta_solver1.solve()
    	phi_solver1.solve()
	q_solver1.solve()

    	eta0.assign(eta1)
    	phi0.assign(phi1)
	q0.assign(q1)

	return (phi0, eta0, q0);

##################################################

