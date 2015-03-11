
from firedrake import *


##################################################

# Intial conditions from exact solutions
def exact_ICs(eta0, phi0, x, kx, ky, Lx, Ly, Ampl, mu, epsilon, speed, is_linear):

	if is_linear == 'True':		# Linear exact solution
		eta0.interpolate(Expression("Ampl*cos(kx*pi*x[0]/Lx)*cos(ky*pi*x[1]/Ly)", kx=kx, ky=ky, Lx=Lx, Ly=Ly, Ampl=Ampl))
    
	elif is_linear == 'False':	# Soliton wave solution
		eta0.interpolate(Expression("1/3.0*c*pow(cosh(0.5*sqrt(c*epsilon/mu)*(x[1]-y0)),-2)", c=speed, mu=mu, epsilon=epsilon, y0=0.5*Ly))
		phi0.interpolate(Expression("2/3.0*sqrt(c*mu/epsilon)*(tanh(0.5*sqrt(c*epsilon/mu)*(x[1]-y0))+1)", c=speed, mu=mu, epsilon=epsilon, y0=0.5*Ly))

	return (eta0, phi0);

##################################################

# Define gravitational potential etaR
def def_etaR(x, H0s, mu, epsilon, hh0, hh1, LL):

	h0 = (hh0-H0s)/(epsilon*H0s)
	h1 = (hh1-H0s)/(epsilon*H0s)
	y1 = LL/(H0s/sqrt(mu))
	y2 = (LL+0.2)/(H0s/sqrt(mu))

	etaR_expr = Expression("(h1-H0)*(0.5*(1+copysign(1.0,y1-x[1])) + 0.25*(1+copysign(1.0,y2-x[1]))*(1+copysign(1.0,x[1]-y1))*(1-(x[1]-y1)/(y2-y1)))", y1=y1, y2=y2, h1=h1, H0=h0)

	return (etaR_expr, h0, h1);

##################################################

# Solve weak form of linear problem eta = etaR for initial eta0
def eta0_eq_etaR(eta,eta0,etaR,gamma):

	aetaR = gamma*eta*dx
	LetaR = gamma*etaR*dx
	etaR_problem = LinearVariationalProblem(aetaR,LetaR,eta0)
	etaR_solver = LinearVariationalSolver(etaR_problem)
	etaR_solver.solve()

	return eta0;

##################################################

