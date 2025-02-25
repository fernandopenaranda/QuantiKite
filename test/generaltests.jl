using Quantica


## 1) GRAPHENE MONOLAYER


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
sdos = dos(c, num_points = 10, num_moments = 10, num_random = 10, num_disorder = 1)
h5gen(h, c, sdos, filename = "testfile.h5")
# Optical conductivity struct
soptcond = conductivity_optical(c,  num_points = 10,  num_moments =10,  num_random = 1,  num_disorder = 1, direction = "yy", temperature = 0.0)
h5gen(h, c, soptcond, filename = "testfile.h5")

## 2) KANE-MELE MODEL
"""
(ref:https://github.com/pablosanjose/Quantica.jl/blob/master/docs/src/tutorial/hamiltonians.md)
 same configuration files `c` as in (1)
"""
function kane_mele_ham()
    SOC(dr) = 0.05 * ifelse(iseven(round(Int, atan(dr[2], dr[1])/(pi/3))), im, -im)
    model = hopping(1, range = neighbors(1)) + 
            hopping((r, dr) ->  SOC(dr); sublats = :A => :A, range = neighbors(2)) +
            hopping((r, dr) -> -SOC(dr); sublats = :B => :B, range = neighbors(2))
    return LatticePresets.honeycomb() |> model
end

h5gen(kane_mele_ham(), c, sdos, filename = "testfile.h5")
h5gen(kane_mele_ham(), c, soptcond, filename = "testfile.h5")


