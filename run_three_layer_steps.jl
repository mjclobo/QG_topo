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



    function topo_rand(h_rms, kt, Lx, Nx)

        # Wavenumber grid
        nkr = Int(Nx / 2 + 1)
        nl = Nx
    
        dk = 2 * pi / Lx
        dl = dk
        
        k = reshape( rfftfreq(Nx, dk * Nx), (nkr, 1) )
        l = reshape( fftfreq(Nx, dl * Nx), (1, nl) )

        k = @. sqrt(k^2 + l^2)

        # Isotropic Gaussian in wavenumber space about mean, kt, with standard deviation, sigma
        # with random Fourier phases
        sigma = sqrt(2) * dk

        seed!(1234)
        hh = exp.(-(k .- kt*(2*pi/Lx)).^2 ./ (2 * sigma^2)) .* exp.(2 * pi * im .* rand(nkr, nl))

        # Recover h from hh
        h = irfft(hh, Nx)

        c = h_rms / sqrt.(mean(h.^2))
        h = c .* h

        return h
    end

    # setting topography
    function topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,type)
        Nx = length(grid_topo.x); Ny = length(grid_topo.y)
        eta_out = zeros(Nx,Ny)
        x = collect(grid_topo.x); y = collect(grid_topo.y)
        for i=1:Nx
            for j=1:Ny
                if type=="eggshell"
                    eta_out[i,j] = (f0/H[end]) * h0 * cos(2*pi*kt*x[i]/Lx) * cos(2*pi*kt*y[j]/Ly)
                elseif type=="sinusoid"
                    eta_out[i,j] = (f0/H[end]) * h0 * cos(2*pi*kt*x[i]/Lx)
                elseif type=="y_slope"
                    eta_out[i,j] = (f0/H[end]) * ((h0*Lx) * ((j-Ny/2)/Ny))
               end
            end
        end

        if type=="rand"
            eta_out = (f0/H[end]) * topo_rand(h0,kt,Lx,Nx)
            # T = eltype(grid_topo)
            # nx, ny = grid_topo.nx, grid_topo.ny
            # h, hx, hy = zeros(T, nx, ny), zeros(T, nx, ny), zeros(T, nx, ny)
            # mfile = matopen("hrand256Km2tk10filtnx32.mat")
            # h = read(mfile, "h")
            # close(mfile)
            # @. h = h*h0*0.5 # ht is rms in random topography. The factor 1/2 is the rms of sin(x)sin(y).
            # hh = rfft(h)
            # hxh = @. im * grid_topo.kr * hh
            # hyh = @. im * grid_topo.l * hh
            # hx = irfft(hxh, nx)
            # hy = irfft(hyh, nx)
            # eta_out = (f0/H[end]) * copy(h)
        end

        return eta_out
    end

    aliased_fraction=1/3; T=Float64;
    grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction, T)
    if topo_type=="eggshell"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"eggshell")
    elseif topo_type=="sinusoid"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sinusoid")
    elseif topo_type=="y_slope"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"y_slope")
    elseif topo_type=="rand"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"rand")
    else
        eta = 0.
    end

    # define the model problem
    prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, g, H, ρ, U, eta=eta,
    μ, β, dt, stepper, linear, aliased_fraction=1/3)

    sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid
    x, y = grid.x, grid.y

    bu_lw = [sqrt(params.g′[1]*(0.5*(H[1]+H[2]))) sqrt(params.g′[2]*(0.5*(H[2]+H[3])))]/f0

    rd_sum = sum(bu_lw)/pi  # where does pi factor come from?

    rd_bulk = sqrt(g*sum(H)*(rho[3]-rho[1])/rho[1])/f0

    # setting initial conditions; does it matter where is 0?
    seed!(1234) # reset of the random number generator for reproducibility
    q₀  = q0_mag * device_array(dev)(randn((grid.nx, grid.ny, nlayers)))
    q₀h = prob.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
    q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

    # q₀[:,:,3] .= 0.0

    MultiLayerQG.set_q!(prob, q₀)

    # output dirs
    filepath = "."

    if topo_type=="y_slope"
        plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/main/"
        plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi/"
        plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
    elseif linear
        plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/main/"
        plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi/"
        plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
    else
        plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/main/"
        plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi/"
        plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi_vert/"
    end
    plotname = "snapshots"
    filename = joinpath(filepath, "2layer.jld2")

    include("./plotting_functions.jl")

    # file management
    if isfile(filename); rm(filename); end
    if !isdir(plotpath_main); mkpath(plotpath_main); end
    if !isdir(plotpath_psi); mkpath(plotpath_psi); end
    if !isdir(plotpath_psi_vert); mkpath(plotpath_psi_vert); end

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
NL1  = [[specE[4][:,:,1]]]
NL2  = [[specE[4][:,:,2]]] 
NL3  = [[specE[4][:,:,3]]] 

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
    eve1,eva1,max_eve1,max_eve_phase1,max_eva1,k_x,k_y,qx1,qy1,rd1 = LinStab.lin_stab(U,V,H,beta,eta,Nx,Ny,rho,f0,g,Float64(Lx),Float64(Ly))
    sigma_LS_all = Array(imag(eva1))
    sigma_LS_mid = sigma_LS_all[:,round(Int,Nx/2)]

    plot_Qy(H,qy1',plotpath_main)

    global bulk_Bu = (real(rd1[2])/Lx)^2
end

# trying to run the model now
startwalltime = time()

cyc = 0
global ell, j = 1, 0

while j<nsteps
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
        push!(NL1,[specE[4][:,:,1]])
        push!(NL2,[specE[4][:,:,2]])
        push!(NL3,[specE[4][:,:,3]])

        # finding vertical structure of instability
        psi_vert1 = abs.(rfft(psi1[:,128]))
        psi_vert2 = abs.(rfft(psi2[:,128]))
        psi_vert3 = abs.(rfft(psi3[:,128]))
        psi_vert = [maximum(psi_vert1)
                    maximum(psi_vert2)
                    maximum(psi_vert3)]
        
        psi_vert = psi_vert./maximum(psi_vert)
    
        # plotting stuff
        global plot_model
        if plot_model==true
            plot_three_layer(tiempo,[KE1 KE2 KE3],[CV32[ell+1] CV52[ell+1] CL1[ell+1] CT[ell+1] NL1[ell+1] NL2[ell+1] NL3[ell+1]],vars.q,grid,kt,h0,plotpath_main,plotname,ell) 

            # Should I also plot the LS most unstable vert structure and the model output for most unstable wavenumber?
            plot_unstable_vert(H,max_eve1,psi_vert,plotpath_psi,plotname,ell)

            plot_layerwise_spectra(grid.kr*Lx/(2*pi),[psi_vert1 psi_vert2 psi_vert3],plotpath_psi_vert,plotname,ell)
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

end

# Get growth rate from exponential fit to upper-layer KE time series.
# if calc_growth_rate==true
#     using Peaks
#     sigma_emp = LinStab.calc_growth(tiempo, [KE1 KE2 KE3 PE32 PE52])
#     a = findmax(CV32[end][1])
#     k_emp = grid.kr[a[2][1]]
# end

if calc_growth_rate==true && perform_ls==true
    plot_growth_rate(k_x[:],sigma_LS_mid,k_emp,sigma_emp,Lx,plotpath_main)
elseif perform_ls==true
    plot_growth_rate(k_x[:],sigma_LS_mid,[],[],Lx,plotpath_main)
end

# save model results, if asked

if save_output
    using CSV

    println("Saving output data to CSV")

    csv_name = "./data/threelayer_" * type*"_gamma"*string(gamma)*"_topo"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) * ".csv"
    # ψ₁, ψ₂ = vars.ψ[:, :, 1], vars.ψ[:, :, 2]

    # should I add streamfunction or PV here?? How would I use them?
    csv_data = Dict("t" => tiempo, "CV32" => CV32, "CV52" => CV52, "KE1" => KE1, "KE2" => KE2, "KE3" => KE3, "Nz" => nlayers, "L" => L,
                    "dt" => dt, "F_profile" => params.F, "beta" => β, "h0" => h0, "kt" => kt, "rho" => rho,
                    "k_ls" => k_x[:], "sigma_ls" => sigma_LS_all,"cfl_set" => cfl_glob)

    CSV.write(csv_name, csv_data)

end
