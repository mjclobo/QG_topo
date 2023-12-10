# This is the file from which the model is ran.
# This is where we set parameters that change.

# whether or not to perform the linear stability analysis
global perform_ls = true

# whether or not to save model output
global save_output = true

# whether or not to plot model output at nsubs timesteps
using PyPlot
global plot_model = true; pygui(false)

# whether or not to calculate growth rate from model output
global calc_growth_rate = true

# type w/ current options: "idealized_Charney", "idealized_Eady", or "real_param_space"
global type = "idealized_Eady"

# set topography w/ current options: "eggshell", "sinusoid"
global topo_type = "eggshell" # "sinusoid" # 

# Hovmoller of streamfunction; also calculates cr and cr_Dopp
global psi_hovm = true

# Magnitude of initial, random QG PV; How to set this??
global q0_mag = 1e-7

# linear or nonlinear model (default: false)
global linear = false

# setting type of run
global run_type = "yrs_to_ss"

if run_type=="power_iter"
    # renormalization parameters
    Rthresh = 0.01        # Threshold energy ratio for renormalization.
    cycles  = 1           # Total number of renormalization cycles.
elseif run_type=="nsteps"
    # number of steps
    nsteps = 30000

    calc_growth_rate = false # for plotting
else
    println("You must choose a valid run type.")
end

# controls ratio of interface densities
gammas = [0.5] # collect(range(0.1,3,5))

# controls ratio of interface shears
alphas = [4.0] # collect(range(1,5,20))

# topo parameters
h0s = [0.] # collect(range(0.,500.,6))      # dimensional topo height 
kts = [50.] # collect(range(1.,50.,6))         # topo wavenumber (no factor of 2pi)

include("./params_three_layer.jl")

# running model
if run_type=="power_iter"
    include("./run_three_layer_power_iter.jl")
elseif run_type=="nsteps"
    include("./run_three_layer_steps.jl")
elseif run_type=="yrs_to_ss"
    include("./run_three_layer_yrs_ss.jl")
else
    println("You must choose a valid run type.")
end

# # Prandlt e-folding scale
# H_t = f0^2 * (Lx/kt)^2 / (g*(rho[3]-rho[2])/rho[1])

# cr_crit = 0.5 * (cr_dopp[1] + cr_dopp[2])

# U32_int = U[1] - (H[1]/2) * ((U[1]-U[2])/((H[1]+H[2])/2))

# H_qy_is_0 = - H[1]/2 - qy1[1] * (((H[1]+H[2])/2))/(qy1[1]-qy1[2])

# U32_eq_cr = U[1] + (H[1]/2 + H_qy_is_0) * ((U[1]-U[2])/((H[1]+H[2])/2))

# cr_int = cr[1] - (H[1]/2) * ((cr[1]-cr[2])/((H[1]+H[2])/2))

# H_cr_is_0 = - H[1]/2 - cr[1] * (((H[1]+H[2])/2))/(cr[1]-cr[2])

# cr_qy_is_0 = cr[1] + (H[1]/2 + H_qy_is_0) * ((cr[1]-cr[2])/((H[1]+H[2])/2))


