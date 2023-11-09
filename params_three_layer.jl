# Parameters for all three-layer model runs;
# both constants and those set by inputs.

using FourierFlows: CPU, TwoDGrid
using Random: seed!, rand
using FFTW

## build basic model

dev = CPU()     # device (CPU)

# numerical params
Ny = Nx = n = 64                # 2D resolution = n²
stepper = "FilteredRK4"         # time stepping scheme
nsubs   = 100                   # number of time-steps for plotting

# physical params
L = Lx = Ly = 400000                    # domain size
μ = 1e-7                                # bottom drag
beta = β = 0                            # the y-gradient of planetary PV

nlayers = 3                 # number of layers
f0 = f₀= 1.24e-4            # Coriolis param [s^-1] 
g = 9.81                    # gravity (in John Mayer voice)
H = [200., 500., 2000.]     # the rest depths of each layer

U = [0.04,0.01,0.0]

# setting density profile as function of gamma
rho = ρ = [0.0, 1025.0, 1025.25]         # the density of each layer
rho1 = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma
rho[1] = ρ[1] = rho1

println("With gamma="*string(gamma)*", your upper layer density is ρ₁="*string(rho1))

# always zero for now, model can't even take this in as a parameter...
V = zeros(nlayers)

# setting cfl and dt; should we switch to a constant dt?...hypothetically this
# should already be constant since U[1] is fixed.
cfl_glob = 0.1
dx = L/Nx
dt = dx*cfl_glob/U[1]

# setting topography
function topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H)
    Nx = length(grid_topo.x); Ny = length(grid_topo.y)
    eta_out = zeros(Nx,Ny)
    x = collect(grid_topo.x); y = collect(grid_topo.y)
    for i=1:Nx
      for j=1:Ny
        eta_out[i,j] =  (f0/H[end]) * h0 * cos(2*pi*kt*x[i]/Lx) * cos(2*pi*kt*y[j]/Ly)
      end
    end
    return eta_out
  end

aliased_fraction=1/3; T=Float64;
grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction, T)
if topo_type=="eggshell"
    eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H)
else
    eta = 0.
end



