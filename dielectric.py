import numpy
import scipy
import sys
import dielectric_tools


if len( sys.argv ) != 5:

    print "usage: python %s h Nx Ny Nz" % sys.argv[0]
    exit(1)

lbkt = 1.0

h  = int( sys.argv[1] )
Nx = int( sys.argv[2] )
Ny = int( sys.argv[3] )
Nz = int( sys.argv[4] )

charge_density = numpy.zeros( [Nz, Ny, Nx], dtype=float )
charge_density_fft_gf = numpy.zeros( [Nz, Ny, Nx/2+1], dtype=complex )
boundary_positions  = numpy.zeros( [0,3], dtype=int )

for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx):
            
            r_sq = (x-Nx/2.)**2 + (y-Ny/2.)**2 + (z-Nz/2.)**2
            
            if r_sq >= 4.**2 and r_sq < 5.**2:
                boundary_positions.resize( [ boundary_positions.shape[0]+1, 3 ] )
                boundary_positions[-1] = [z, y, x]
                
            if r_sq < 5.**2:
                charge_density[z,y,x] = 1.
#charge_density[0,0,0] = 1.0

node_boundary = numpy.zeros( [Nz, Ny, Nx], dtype=int )

for i in boundary_positions:
    node_boundary[i[0], i[1], i[2]] = 1
    
dielectric_tools.write_scalar_vtk( charge_density, "charge_density.vtk" )
dielectric_tools.write_scalar_vtk( node_boundary, "boundary.vtk" )

charge_density_fft = numpy.fft.rfftn(charge_density)

dielectric_tools.write_complex_vtk( charge_density_fft, "charge_density_fft.vtk" )

for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx/2+1):
            if (x,y,z) == (0,0,0):
                charge_density_fft_gf[z,y,x] = 0.0
            else:
                charge_density_fft_gf[z,y,x] = charge_density_fft[z,y,x] * (-4.0) * numpy.pi * lbkt * h**2 * 0.5 / ( numpy.cos( 2.0 * numpy.pi * x / Nx ) + numpy.cos( 2.0 * numpy.pi * y / Ny ) + numpy.cos( 2.0 * numpy.pi * z / Nz ) - 3.0 )

dielectric_tools.write_complex_vtk( charge_density_fft_gf, "charge_density_fft_gf.vtk" )

charge_density_fft_gf_ifft = numpy.fft.irfftn( charge_density_fft_gf )

dielectric_tools.write_scalar_vtk( charge_density_fft_gf_ifft, "charge_density_fft_gf_ifft.vtk" )
