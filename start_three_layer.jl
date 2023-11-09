# This is the file from which the model is ran.
# This is where we set parameters that change.

# controls ratios of two upper layer densities spanning either side of ``critical'' density (see equation)
gamma = 1.5

# whether or not to perform the linear stability analysis
perform_ls = true

# whether or not to save model output
save_output = true

# whether or not to plot model output at nsubs timesteps
global plot_model = true

# whether or not to calculate growth rate from model output
calc_growth_rate = true

# type w/ current options: "idealized_Charney", "idealized_Eady", or "real_param_space"
type = "idealized_Eady"

# set topography w/ current option(s): "eggshell"
# currently the height is fixed, but maybe it's worth varying this height for
# a given Charney-type profile in the future..
topo_type = "eggshell"

# Magnitude of initial, random QG PV; How to set this??
q0_mag = 1e-8

# topo
h0 = 0.       # dimensional topo height; constant for now
kt = 6.         # topo wavenumber (no factor of 2pi); constant for now

# renormalization parameters
Rthresh = 0.001        # Threshold energy ratio for renormalization.
cycles  = 5           # Total number of renormalization cycles.

# loading shared params
include("./params_three_layer.jl")

# running model
include("./run_three_layer.jl")




