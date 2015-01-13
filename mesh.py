
from firedrake import *


##################################################

def create_mesh(Nx, Ny, Lx, Ly):

	m = UnitIntervalMesh(Nx)            # For quadrilateral mesh
	mesh = ExtrudedMesh(m, layers=Ny)
	#mesh = UnitSquareMesh(Nx,Ny)       # For triangular mesh
	coords = mesh.coordinates
	coords.dat.data[:,0] = Lx*coords.dat.data[:,0]
	coords.dat.data[:,1] = Ly*coords.dat.data[:,1]

	return (mesh, coords);

##################################################

def stretch_mesh(mesh, coords, Nx, Ny, Lx, d):

	L = sqrt(d**2 - (0.5*Lx)**2)
	slope = L/(0.5*Lx)
	dy = [0]*(Nx+1)*(Ny+1)

	for i in range(len(coords.dat.data)):
		if coords.dat.data[i,0]<0.5*Lx:
			dy[i] = slope*(coords.dat.data[i,0]-0.5*Lx)/Ny
		else:
		        dy[i] = slope*(0.5*Lx-coords.dat.data[i,0])/Ny
		coords.dat.data[i,1] = coords.dat.data[i,1] + i%(Ny+1)*dy[i]

	return coords;

##################################################

