import numpy


def write_scalar_vtk( X, filename, name="python_data", origin=(0., 0., 0.), spacing=(1., 1., 1.) ):

    dim = numpy.ones( 3, dtype=int )
    
    for i in numpy.arange( 0, X.ndim ):
        dim[i] = X.shape[i]

    vtkout = open( filename, 'w' )
    
    vtkout.write( "# vtk DataFile Version 2.0\n" )
    vtkout.write( "%s\n" % ( name, ) )
    vtkout.write( "ASCII\n\n" )
    
    vtkout.write( "DATASET STRUCTURED_POINTS\n" )
    vtkout.write( "DIMENSIONS %s %s %s\n" % ( str(dim[2]),     str(dim[1]),     str(dim[0])     ) )
    vtkout.write( "ORIGIN %s %s %s\n"     % ( str(origin[0]),  str(origin[1]),   str(origin[2])  ) )
    vtkout.write( "SPACING %s %s %s\n\n"  % ( str(spacing[0]), str(spacing[1]), str(spacing[2]) ) )
    
    vtkout.write( "POINT_DATA %s\n"       % ( str( dim[0]*dim[1]*dim[2] ) ) )
    vtkout.write( "SCALARS %s FLOAT 1\n"  % ( name, ) )
    vtkout.write( "LOOKUP_TABLE default\n" )
    
    for z in numpy.arange(0, dim[0]):
        for y in numpy.arange(0, dim[1]):
            for x in numpy.arange(0, dim[2]):
                vtkout.write( "%s " % ( str(X[z,y,x]), ) )
    
    vtkout.close()


def write_complex_vtk( X, filename, name="python_data", origin=(0., 0., 0.), spacing=(1., 1., 1.) ):

    dim = numpy.ones( 3, dtype=int )
    
    for i in numpy.arange( 0, X.ndim ):
        dim[i] = X.shape[i]

    vtkout = open( filename, 'w' )
    
    vtkout.write( "# vtk DataFile Version 2.0\n" )
    vtkout.write( "%s\n" % ( name, ) )
    vtkout.write( "ASCII\n\n" )
    
    vtkout.write( "DATASET STRUCTURED_POINTS\n" )
    vtkout.write( "DIMENSIONS %s %s %s\n" % ( str(dim[2]),     str(dim[1]),     str(dim[0])     ) )
    vtkout.write( "ORIGIN %s %s %s\n"     % ( str(origin[0]),  str(origin[1]),   str(origin[2])  ) )
    vtkout.write( "SPACING %s %s %s\n\n"  % ( str(spacing[0]), str(spacing[1]), str(spacing[2]) ) )
    
    vtkout.write( "POINT_DATA %s\n"       % ( str( dim[0]*dim[1]*dim[2] ) ) )
    vtkout.write( "SCALARS %s FLOAT 2\n"  % ( name, ) )
    vtkout.write( "LOOKUP_TABLE default\n" )
    
    for z in numpy.arange(0, dim[0]):
        for y in numpy.arange(0, dim[1]):
            for x in numpy.arange(0, dim[2]):
                vtkout.write( "%s %s " % ( str(X[z,y,x].imag), str(X[z,y,x].real) ) )
    
    vtkout.close()
