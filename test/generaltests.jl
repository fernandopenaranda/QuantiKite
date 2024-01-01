using Quantica

h = Quantica.HamiltonianPresets.graphene()
# Configuration
c = configuration(
    divisions = (1, 1), 
    unitcells =(32, 32), 
    boundaries = ("periodic", "periodic"), 
    precision = 2, 
    spectrum_range = nothing,
    angles = (0., 0.)
    )
# dos 
s = dos(c, num_points = 10, num_moments = 10, num_random = 10, num_disorder = 1)
h5gen(h, c, s, filename = "testfile.h5")
# Optical conductivity struct
s = conductivity_optical(c,  num_points = 10,  num_moments =10,  num_random = 1,  num_disorder = 1, direction = "yy", temperature = 0.0)
h5gen(h, c, s, filename = "testfile.h5")
