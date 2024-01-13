# Parameters for all three-layer model runs;
# both constants and those set by inputs.

using FourierFlows: CPU, TwoDGrid
using Random: seed!, rand
using FFTW
using GeophysicalFlows, Printf, FFTW, LinearAlgebra, Statistics, LaTeXStrings, CSV, Peaks
using Random: seed!

# using Pkg
# Pkg.develop(path="../gfjl/GeophysicalFlows.jl")


include("../LinStab/mjcl_stab.jl")
using .LinStab

## build basic model

dev = CPU()     # device (CPU)

# numerical params
Ny = Nx = n = 64                # 2D resolution = n²
stepper = "FilteredRK4"         # time stepping scheme
nsubs   = 100                   # number of time-steps for plotting; for nsteps this is set in run_three_layer_nsteps.jl!!!

# physical params
L = Lx = Ly = 1000.e3                   # domain size [m]
beta = β = 0 # 1.9e-11 # 1.14052e-11              # the y-gradient of planetary PV

nlayers = 3                 # number of layers
f0 = f₀= 8.3e-5  # 1.0e-4            # Coriolis param [s^-1] 
g = 9.81                    # gravity
H = [500., 1000., 2500.]     # the rest depths of each layer; [250., 500., 3000.] 

H32 = 0.5 * (H[1] + H[2])
H52 = 0.5 * (H[2] + H[3])

U = [0.05,0.0255,0.0]

# setting base density profile
rho = ρ = [0.0, 1025.0, 1025.75]         # the density of each layer

# rho[3] = rho[2] + rho[2] - rho[1]

# println("With gamma="*string(gamma)*", your upper layer density is ρ₁="*string(rho1))

# from Wenda/s code
# rho = ρ = [1026.50720363, 1027.50749062, 1027.91389085]
# H = [440., 1760., 1760.]
# U = [0.02659368, 0.00457271, 0.00013535]

# gamma = 2.5
# alpha = 0.5

# U[2] = U[1] - alpha*(U[1]-U[3])

# rho[1] = ρ[1] = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma

# always zero for now, model can't even take this in as a parameter...
V = zeros(nlayers)

# Linear bottom drag
# μ = 1e-7                                # bottom drag
μ =  (86400*20)^-1           # a typical value of 5-10 m in the ocean [Weatherly and Martin, 1978; Perlin et al., 2007]

# setting cfl and dt; should we switch to a constant dt?...hypothetically this
# should already be constant since U[1] is fixed.
cfl_glob = 0.0125 # 0.1 is standard so far; 0.5 for 128
dx = L/Nx
dt = dx*cfl_glob/U[1]     # /5

## alternate way to set density values
# N2 = [2*10^-6 2*10^-5 2*10^-6]

# rho0 = 1025.
# rho12 = rho0
# rho32 = rho12 + (rho0*H[1]*N2[1])/g
# rho52 = rho32 + (rho0*H[2]*N2[2])/g
# rho72 = rho52 + (rho0*H[3]*N2[3])/g

# rho = [0.5*(rho12+rho32), 0.5*(rho32+rho52), 0.5*(rho52+rho72)]

