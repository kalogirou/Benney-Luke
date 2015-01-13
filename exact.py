
from firedrake import *


##################################################

def exact_solutions(x, t, kx, ky, k2, Lx, Ly, omega, Ampl, mu, epsilon, speed, is_linear):

	if is_linear == 'True':
		# Exact linear solution
		expr_phi = Expression("-Ampl/w/(1+0.5*mu*k2)*sin(w*t)*cos(kx*pi*x[0]/Lx)*cos(ky*pi*x[1]/Ly)", t=t, kx=kx, ky=ky, Lx=Lx, Ly=Ly, w=omega, Ampl=Ampl, mu=mu, k2=k2)
		expr_eta = Expression("Ampl*cos(w*t)*cos(kx*pi*x[0]/Lx)*cos(ky*pi*x[1]/Ly)", t=t, kx=kx, ky=ky, Lx=Lx, Ly=Ly, w=omega, Ampl=Ampl)
    
	elif is_linear == 'False':
		# Exact soliton solution
		expr_eta = Expression("1/3.0*c*pow(cosh(0.5*sqrt(c*epsilon/mu)*(x[1]-y0-t-epsilon*c*t/6.0)),-2)", t=t, c=speed, epsilon=epsilon, mu=mu, y0=0.5*Ly)
		expr_phi = Expression("2/3.0*sqrt(c*mu/epsilon)*(tanh(0.5*sqrt(c*epsilon/mu)*(x[1]-y0-t-epsilon*c*t/6.0))+1)", t=t, c=speed, epsilon=epsilon, mu=mu, y0=0.5*Ly)

	return (expr_phi,expr_eta);

##################################################

def update_exact(ex_phi, ex_eta, expr_phi, expr_eta, t, Ly, epsilon, c, is_linear):

	if is_linear == 'True':
		expr_eta.t = t					# Update exact solutions in time
    		expr_phi.t = t
    		
	elif is_linear == 'False':				# Soliton reflection on solid wall
		if t>1.5*Ly/(1+epsilon*c):
			expr_eta.t = t-1.5*Ly/(1+epsilon*c)
		        expr_phi.t = t-1.5*Ly/(1+epsilon*c)
			expr_eta.y0 = 0
			expr_phi.y0 = 0
	    	elif t>0.5*Ly/(1+epsilon*c):
			expr_eta.t = 0.5*Ly/(1+epsilon*c)-t
        		expr_phi.t = 0.5*Ly/(1+epsilon*c)-t
			expr_eta.y0 = Ly
			expr_phi.y0 = Ly
    		else:
			expr_eta.t = t
    			expr_phi.t = t

    	ex_phi.interpolate(expr_phi)
    	ex_eta.interpolate(expr_eta)

	return (ex_phi, ex_eta);

##################################################

