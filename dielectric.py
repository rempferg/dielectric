import numpy
import scipy
import sys
import dielectric_tools


if len( sys.argv ) != 5:

  print "usage: python %s h Nx Ny Nz" % sys.argv[0]
  exit(1)
  
  
h  = int( sys.argv[1] )
Nx = int( sys.argv[2] )
Ny = int( sys.argv[3] )
Nz = int( sys.argv[4] )

charge_density = numpy.zeros( [Nz, Ny, Nx], dtype=float )
boundary_positions  = numpy.zeros( [0,3], dtype=int )

for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx):
            
            r_sq = (x-Nx/2.)**2 + (y-Ny/2.)**2 + (z-Nz/2.)**2
            
            if r_sq >= 4.**2 and r_sq < 5.**2:
                boundary_positions.resize( [ boundary_positions.shape[0]+1, 3 ] )
                boundary_positions[-1] = [z, y, x]

node_boundary = numpy.zeros( [Nz, Ny, Nx], dtype=int )

for i in boundary_positions:
    node_boundary[i[0], i[1], i[2]] = 1

dielectric_tools.write_scalar_vtk( node_boundary, "charge_density.vtk" )

charge_density_fft = numpy.fft.rfftn(charge_density)

dielectric_tools.write_scalar_vtk( node_boundary, "charge_density_fft.vtk" )
