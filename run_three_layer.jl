# This is the kernel for running QG models using GeophysicalFlows.jl.
# You pass in all relevant params, and this script executes n_iter
# iterations of the power renormalization method.
# Still an open question: How do we know what initial conditions are
# appropriate, given that the condition for reiteration is 100x I.C.?
# Maybe it is worth running different orders of magnitude of I.C.s and
# seeing if/how this affects the instability results.

## load packages

using GeophysicalFlows, Plots, Printf, FFTW, LinearAlgebra, Statistics, LaTeXStrings

using Random: seed!

# define the model problem
prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, g, H, ρ, U, eta=eta,
μ, β, dt, stepper, aliased_fraction=0)

sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid
x, y = grid.x, grid.y

# setting initial conditions; does it matter where is 0?
seed!(1234) # reset of the random number generator for reproducibility
q₀  = q0_mag * device_array(dev)(randn((grid.nx, grid.ny, nlayers)))
q₀h = prob.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

MultiLayerQG.set_q!(prob, q₀)

# output dirs
filepath = "."
plotpath_main = "./figs/plots_3layer"*"_gamma"*string(gamma)*"_h0"*string(Int(h0))*"_res" * string(Int(Nx)) *"/main/"
plotpath_psi  = "./figs/plots_3layer"*"_gamma"*string(gamma)*"_h0"*string(Int(h0))*"_res" * string(Int(Nx)) *"/psi/"
plotname = "snapshots"
filename = joinpath(filepath, "2layer.jld2")

include("./plotting_functions.jl")

# file management
if isfile(filename); rm(filename); end
if !isdir(plotpath_main); mkpath(plotpath_main); end
if !isdir(plotpath_psi); mkpath(plotpath_psi); end

# # ``create output'' (?)
# get_sol(prob) = prob.sol # extracts the Fourier-transformed solution

# function get_u(prob)
#     sol, params, vars, grid = prob.sol, prob.params, prob.vars, prob.grid

#     @. vars.qh = sol
#     streamfunctionfrompv!(vars.ψh, vars.qh, params, grid)
#     @. vars.uh = -im * grid.l * vars.ψh
#     invtransform!(vars.u, vars.uh, params)

#     return vars.u
# end

# out = Output(prob, filename, (:sol, get_sol), (:u, get_u))

# savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), 0)

# defining initial diagnostics
E = MultiLayerQG.energies(prob)
specE = MultiLayerQG.spectralfluxes(prob)

# variables for plotting, that will be pushed to
tiempo = [0.]
KE1 = [E[1][1]]/H[1]
KE2 = [E[1][2]]/H[2]
KE3 = [E[1][3]]/H[3]

PE32 = [E[2][1]]/((H[1]+H[2])/2)
PE52 = [E[2][2]]/((H[2]+H[3])/2)

CV32 = [[specE[2][:,:,1]]]
CV52 = [[specE[2][:,:,2]]]
CL1  = [[specE[1][:,1]]]
CT   = [[specE[3]]]

# psi1 = Array(vars.ψ[:, :, 1])
# psi2 = Array(vars.ψ[:, :, 2])
# psi3 = Array(vars.ψ[:, :, 3])

# initial KE of upper layer, for renormalization
KE1_0 = E[1][1][1]

# Perform linear stability analysis, if asked
if perform_ls==true
    # load module
    include("../LinStab/mjcl_stab.jl")
    using .LinStab

    # perform stability analysis
    eta=0;
    eve1,eva1,max_eve1,max_eva1,k_x,k_y,qx1,qy1,rd1 = LinStab.lin_stab(U,V,H,beta,eta,Nx,Ny,rho,f0,g,Float64(Lx),Float64(Ly))
    sigma_LS_all = Array(imag(eva1))
    sigma_LS_mid = sigma_LS_all[:,round(Int,Nx/2)]

    plot_Qy(H,qy1',plotpath_main)

end

# trying to run the model now
startwalltime = time()

cyc = 0
global ell, j = 1, 0

while cyc<cycles
    global ell, j 

    ##########################
    stepforward!(prob)
    MultiLayerQG.updatevars!(prob)
    local E = MultiLayerQG.energies(prob)
    local specE = MultiLayerQG.spectralfluxes(prob)

    if j % nsubs == 0
        # updating variables to plot
        psi1 = vars.ψ[:, :, 1]
        psi2 = vars.ψ[:, :, 2]
        psi3 = vars.ψ[:, :, 3]
    
        q1 = transpose(vars.q[:, :, 1])
        q2 = transpose(vars.q[:, :, 2])
        q3 = transpose(vars.q[:, :, 3])
    
        # push to time
        push!(tiempo,clock.t)

        # push to layer-wise KE
        push!(KE1,E[1][1]/H[1])
        push!(KE2,E[1][2]/H[2])
        push!(KE3,E[1][3]/H[3])

        # push to interface potential energies
        push!(PE32,E[2][1]/((H[1]+H[2])/2))
        push!(PE52,E[2][2]/((H[2]+H[3])/2))

        # push to spectral flux terms
        push!(CV32,[specE[2][:,:,1]])
        push!(CV52,[specE[2][:,:,2]])
        push!(CL1,[specE[1][:,1]])
        push!(CT,[specE[3]])

        # finding vertical structure of instability
        psi_vert = [maximum(abs.(rfft(psi1[:,32])))
                    maximum(abs.(rfft(psi2[:,32])))
                    maximum(abs.(rfft(psi3[:,32])))]
        
        psi_vert = psi_vert./maximum(psi_vert)
    
        # plotting stuff
        global plot_model
        if plot_model==true
            plot_three_layer(tiempo,[KE1 KE2 KE3],[CV32[ell+1] CV52[ell+1] CL1[ell+1] CT[ell+1]],vars.q,grid,kt,h0,plotpath_main,plotname,ell)

            # Should I also plot the LS most unstable vert structure and the model output for most unstable wavenumber?
            plot_unstable_vert(H,max_eve1,psi_vert,plotpath_psi,plotname,ell)
        end

        # increase counter
        global ell+=1
    
        # reading out stats
        cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
    
        log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1: %.3e, KE2: %.3e, PE: %.3e, walltime: %.2f min",
                        clock.step, clock.t, cfl, E[1][1], E[1][2], E[2][1], (time()-startwalltime)/60)
        
        println(log)

    end

    global j+=1
    ###############################

    R = KE1_0/E[1][1] # renormalize (KE_{ref}/KE).
    if R<Rthresh
        MultiLayerQG.set_q!(prob, vars.q*R)
        global cyc += 1
        println("")
        println("Completed renormalization cycle", " ", cyc, "/", cycles)
        println("")
    end

end

# Get growth rate from exponential fit to upper-layer KE time series.
if calc_growth_rate==true
    using Peaks
    sigma_emp = LinStab.calc_growth(tiempo, [KE1 KE2 KE3 PE32 PE52])
    a = findmax(CV32[end][1])
    k_emp = grid.kr[a[2][1]]
end

if calc_growth_rate==true && perform_ls==true
    plot_growth_rate(k_x[:],sigma_LS_mid,k_emp,sigma_emp,Lx,plotpath_main)
end

# save model results, if asked
if save_output
    using CSV

    println("Saving output data to CSV")

    csv_name = "./data/threelayer_" * type*"_gamma"*string(gamma)*"_topo"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) * ".csv"
    # ψ₁, ψ₂ = vars.ψ[:, :, 1], vars.ψ[:, :, 2]

    # should I add streamfunction or PV here?? How would I use them?
    csv_data = Dict("t" => tiempo, "CV32" => CV32, "CV52" => CV52, "KE1" => KE1, "KE2" => KE2, "KE3" => KE3, "Nz" => nlayers, "L" => L,
                    "dt" => dt, "F_profile" => params.F, "beta" => β, "h0" => h0, "kt" => kt, "sigma_emp" => sigma_emp,
                    "k_emp" => k_emp, "k_ls" => k_x[:], "sigma_ls" => sigma_LS_all,"cfl_set" => cfl_glob)

    CSV.write(csv_name, csv_data)

end
