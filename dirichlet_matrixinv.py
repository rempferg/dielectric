import numpy
import scipy
import sys
import dielectric_tools


#check for proper usage
if len( sys.argv ) != 5:

    print "usage: python %s h Nx Ny Nz" % sys.argv[0]
    exit(1)


#setup parameters
lbkt = 1.0

h  = int( sys.argv[1] )
Nx = int( sys.argv[2] )
Ny = int( sys.argv[3] )
Nz = int( sys.argv[4] )


#allocate memory
charge_potential = numpy.zeros( [Nz, Ny, Nx], dtype=float )
tmp_grid_fft = numpy.zeros( [Nz, Ny, Nx/2+1], dtype=complex )
boundary_positions  = numpy.zeros( [0,3], dtype=int )


#initialize boundary geometry
for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx):
            
            r_sq = (x-Nx/2.)**2 + (y-Ny/2.)**2 + (z-Nz/2.)**2
            
            if r_sq >= 8.**2 and r_sq < 9.**2:
                boundary_positions.resize( [ boundary_positions.shape[0]+1, 3 ] )
                boundary_positions[-1] = [z, y, x]


#allocate boundary potential/charge array and element interaction matrix
boundary_potential_charge = numpy.zeros( len(boundary_positions), dtype=float )
boundary_interaction_matrix = numpy.matrix( numpy.zeros( [ len(boundary_positions), len(boundary_positions) ], dtype=float ), dtype=float, copy=False )


#calculate greensfunction
charge_potential[0,0,0] = 1.0

tmp_grid_fft = numpy.fft.rfftn(charge_potential)

for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx/2+1):
        
            if (x,y,z) == (0,0,0):
                tmp_grid_fft[z,y,x] = 0.0
            else:
                tmp_grid_fft[z,y,x] = tmp_grid_fft[z,y,x] * (-4.0) * numpy.pi * lbkt * h**2 * 0.5 / ( numpy.cos( 2.0 * numpy.pi * x / Nx ) + numpy.cos( 2.0 * numpy.pi * y / Ny ) + numpy.cos( 2.0 * numpy.pi * z / Nz ) - 3.0 )

charge_potential = numpy.fft.irfftn( tmp_grid_fft )


#populate boundary element interaction matrix and invert
for i in numpy.arange( 0, len(boundary_positions) ):
    for k in numpy.arange( 0, len(boundary_positions) ):
    
        d = numpy.subtract( boundary_positions[i], boundary_positions[k] )
        d = numpy.mod( d, [ Nx, Ny, Nz ] )
        
        boundary_interaction_matrix[i,k] = charge_potential[ d[2], d[1], d[0] ]

print "Matrix has size", boundary_interaction_matrix.shape

boundary_interaction_matrix_inv = numpy.linalg.inv( boundary_interaction_matrix )


#calculate potential without boundaries
for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx):
        
            charge_potential[z,y,x] = x


#collect vector containing boundary potential
for i in numpy.arange( 0, len(boundary_positions) ):

    boundary_potential_charge[i] = charge_potential[ boundary_positions[i,2], boundary_positions[i,1], boundary_positions[i,0] ]


#calculate influence charge density
boundary_potential_charge = -boundary_interaction_matrix_inv.dot(boundary_potential_charge)


#expand it into whole field and output
charge_potential = numpy.zeros( [Nz, Ny, Nx], dtype=float )

for i in numpy.arange( 0, len(boundary_positions) ):
    
    charge_potential[ boundary_positions[i,2], boundary_positions[i,1], boundary_positions[i,0] ] = boundary_potential_charge[0,i]

dielectric_tools.write_scalar_vtk( charge_potential, "influence_charge.vtk" )


#calculate potential of influence charge
tmp_grid_fft = numpy.fft.rfftn(charge_potential)

for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx/2+1):
        
            if (x,y,z) == (0,0,0):
                tmp_grid_fft[z,y,x] = 0.0
            else:
                tmp_grid_fft[z,y,x] = tmp_grid_fft[z,y,x] * (-4.0) * numpy.pi * lbkt * h**2 * 0.5 / ( numpy.cos( 2.0 * numpy.pi * x / Nx ) + numpy.cos( 2.0 * numpy.pi * y / Ny ) + numpy.cos( 2.0 * numpy.pi * z / Nz ) - 3.0 )

charge_potential = numpy.fft.irfftn( tmp_grid_fft )


#add potential
for z in numpy.arange(0, Nz):
    for y in numpy.arange(0, Ny):
        for x in numpy.arange(0, Nx):
        
            charge_potential[z,y,x] += x

dielectric_tools.write_scalar_vtk( charge_potential, "potential.vtk" )
