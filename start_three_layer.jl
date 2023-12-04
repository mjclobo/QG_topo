# This is the file from which the model is ran.
# This is where we set parameters that change.

# controls ratios of two upper layer densities spanning either side of ``critical'' density (see equation)
gamma = 1.1      # 5.0 is a good number to use; 1.0 is stable with Callies parameters (but with double upper shear)
alpha = 2.2

# topo
h0 = 50.      # dimensional topo height; constant for now
kt = 6.         # topo wavenumber (no factor of 2pi); constant for now

# whether or not to perform the linear stability analysis
perform_ls = true

# whether or not to save model output
save_output = true

# whether or not to plot model output at nsubs timesteps
using PyPlot
global plot_model = false; pygui(false)

# whether or not to calculate growth rate from model output
calc_growth_rate = true

# type w/ current options: "idealized_Charney", "idealized_Eady", or "real_param_space"
type = "idealized_Eady"

# set topography w/ current option(s): "eggshell"
# currently the height is fixed, but maybe it's worth varying this height for
# a given Charney-type profile in the future..
topo_type = "eggshell" # "sinusoid" # 

# Hovmoller of streamfunction
psi_hovm = true

# Magnitude of initial, random QG PV; How to set this??
q0_mag = 1e-7

# linear or nonlinear model (default: false)
linear = false

# setting type of run
run_type = "power_iter"

if run_type=="power_iter"
    # renormalization parameters
    Rthresh = 0.01        # Threshold energy ratio for renormalization.
    cycles  = 3           # Total number of renormalization cycles.
elseif run_type=="nsteps"
    # number of steps
    nsteps = 30000

    calc_growth_rate = false # for plotting
else
    println("You must choose a valid run type.")
end

# loading shared params
include("./params_three_layer.jl")

# running model
if run_type=="power_iter"
    include("./run_three_layer_power_iter.jl")
elseif run_type=="nsteps"
    include("./run_three_layer_nsteps.jl")
else
    println("You must choose a valid run type.")
end



using PyPlot; matplotlib[:rcParams]["axes.unicode_minus"]=false

fig,ax = PyPlot.subplots(1,3,figsize=(12,8));
fig.tight_layout(pad=5.0)
ax1=ax[1]; ax2=ax[2]; ax3=ax[3];

cr1_dopp = cr_dopp[1]; cr2_dopp = cr_dopp[2]; cr3_dopp = cr_dopp[3];

pc1=ax1.pcolormesh(x/1.e3,t_hovm,psi1_ot'./ maximum(abs.(psi1_ot)',dims=2))
ax1.tick_params(labelsize=10.)
ax1.set_xlabel(L"x [km]", fontsize=18.)
ax1.set_ylabel(L"t [s]", fontsize=18.)
ax1.set_title(L"\psi_1 [\mathrm{norm}]", fontsize = 20.)
ax1.set_ylim([t_hovm[1],t_hovm[end]])
ax1.text(0.2,0.1,L"c_{r,1} = " * string(Int(round(cr1_dopp*1000))/1000), transform=ax1.transAxes,fontsize=16.)
fig.colorbar(pc1)
pc2=ax2.pcolormesh(x/1.e3,t_hovm,psi2_ot'./ maximum(abs.(psi2_ot)',dims=2))
ax2.tick_params(labelsize=10.)
ax2.set_xlabel(L"x [km]", fontsize=18.)
ax2.set_ylabel(L"t [s]", fontsize=18.)
ax2.set_title(L"\psi_2 [\mathrm{norm}]", fontsize = 20.)
ax2.set_ylim([t_hovm[1],t_hovm[end]])
ax2.text(0.2,0.1,L"c_{r,2} = " * string(Int(round(cr2_dopp*1000))/1000), transform=ax2.transAxes,fontsize=16.)
fig.colorbar(pc2)
pc3=ax3.pcolormesh(x/1.e3,t_hovm,psi3_ot'./ maximum(abs.(psi3_ot)',dims=2))
ax3.tick_params(labelsize=10.)
ax3.set_xlabel(L"x [km]", fontsize=18.)
ax3.set_ylabel(L"t [s]", fontsize=18.)
ax3.set_title(L"\psi_3 [\mathrm{norm}]", fontsize = 20.)
ax3.set_ylim([t_hovm[1],t_hovm[end]])
ax3.text(0.2,0.1,L"c_{r,3} = " * string(Int(round(cr3_dopp*1000))/1000), transform=ax3.transAxes,fontsize=16.)
fig.colorbar(pc3)

# PyPlot.savefig(fig,plotpath_main*"../psi_hovm_gamma"*string(gamma)*".png")

# Dimensionless parameters
Ri = [((0.5*(H[1]+H[2])*g*(rho[2]-rho[1])/rho[1])/((U[1]-U[2])^2))
        ((0.5*(H[2]+H[3])*g*(rho[3]-rho[2])/rho[1])/((U[2]-U[3])^2))]

Bu_11 = (H[1]*g*(rho[2]-rho[1])/rho[1])/f0^2 * Lx^-2
Bu_12 = (H[2]*g*(rho[2]-rho[1])/rho[1])/f0^2 * Lx^-2
Bu_22 = (H[2]*g*(rho[3]-rho[2])/rho[1])/f0^2 * Lx^-2
Bu_23 = (H[3]*g*(rho[3]-rho[2])/rho[1])/f0^2 * Lx^-2

Bu = [Bu_11, Bu_12, Bu_22, Bu_23]

inv2_Bu = sqrt.(Bu).^-1

# Prandlt e-folding scale
H_t = f0^2 * (Lx/kt)^2 / (g*(rho[3]-rho[2])/rho[1])

cr_crit = 0.5 * (cr_dopp[1] + cr_dopp[2])

U32_int = U[1] - (H[1]/2) * ((U[1]-U[2])/((H[1]+H[2])/2))

H_qy_is_0 = - H[1]/2 - qy1[1] * (((H[1]+H[2])/2))/(qy1[1]-qy1[2])

U32_eq_cr = U[1] + (H[1]/2 + H_qy_is_0) * ((U[1]-U[2])/((H[1]+H[2])/2))

cr_int = cr[1] - (H[1]/2) * ((cr[1]-cr[2])/((H[1]+H[2])/2))

H_cr_is_0 = - H[1]/2 - cr[1] * (((H[1]+H[2])/2))/(cr[1]-cr[2])

cr_qy_is_0 = cr[1] + (H[1]/2 + H_qy_is_0) * ((cr[1]-cr[2])/((H[1]+H[2])/2))


