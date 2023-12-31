# This is the kernel for running QG models using GeophysicalFlows.jl.
# You pass in all relevant params, and this script executes n_iter
# iterations of the power renormalization method.
# Still an open question: How do we know what initial conditions are
# appropriate, given that the condition for reiteration is 100x I.C.?
# Maybe it is worth running different orders of magnitude of I.C.s and
# seeing if/how this affects the instability results.

## load packages

nsubs = 100

for gamma=gammas; for alpha=alphas; for h0=h0s; for kt=kts

    # change variable params
    U[2] = (U[1] + alpha * U[3] * (H32/H52))/(1 + alpha * (H32/H52))

    # rho1 = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma
    rho1 = rho[2] - gamma * (rho[3] - rho[2]) * (H32/H52)

    rho[1] = ρ[1] = rho1

    # setting topography
    function topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,type)
        Nx = length(grid_topo.x); Ny = length(grid_topo.y)
        eta_out = zeros(Nx,Ny)
        x = collect(grid_topo.x); y = collect(grid_topo.y)
        for i=1:Nx
            for j=1:Ny
                if type=="eggshell"
                    eta_out[i,j] = (f0/sum(H)) * h0 * cos(2*pi*kt*x[i]/Lx) * cos(2*pi*kt*y[j]/Ly)
                elseif type=="sinusoid"
                    eta_out[i,j] = (f0/sum(H)) * h0 * cos(2*pi*kt*x[i]/Lx)
                end
            end
        end
        return eta_out
    end

    aliased_fraction=1/3; T=Float64;
    grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction, T)
    if topo_type=="eggshell"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"eggshell")
    elseif topo_type=="sinusoid"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sinusoid")
    else
        eta = 0.
    end

    nν=4
    ν=3.0e9

    # define the model problem
    prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, g, H, ρ, U, eta=eta, nν, ν,
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

    if linear
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

    # Perform linear stability analysis, if asked
    if perform_ls==true
        # perform stability analysis
        eta=0;
        eve1,eva1,max_eve1,max_eve_phase1,max_eva1,k_x,k_y,qx1,qy1,rd1 = LinStab.lin_stab(U,V,H,beta,eta,Nx,Ny,rho,f0,g,Float64(Lx),Float64(Ly))
        sigma_LS_all = Array(imag(eva1))
        sigma_LS_mid = sigma_LS_all[:,round(Int,Nx/2)]

        # plot_Qy(H,qy1',plotpath_main)

        global bulk_Bu = (real(rd1[2])/Lx)^2
    end

    # trying to run the model now
    startwalltime = time()

    global ell, j = 1, 0
    global t_hovm = nothing
    global yr_cnt = 1
    global ss_yr = false
    global ss_yr_cnt = 0
    global ss_yr_max = 5
    global KE_thresh = 1.05      # if KE at end of year over KE at beginning of yr is less than this, then model is in steady-state

    while ss_yr_cnt < ss_yr_max
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
        
            if isnothing(t_hovm)
                global psi1_ot = psi1[:,Int(round(Nx/2))]
                global psi2_ot = psi2[:,Int(round(Nx/2))]
                global psi3_ot = psi3[:,Int(round(Nx/2))]
                global t_hovm = Array([clock.t])

                global psi_vert1 = abs.(rfft(psi1[:,32]))
                global psi_vert2 = abs.(rfft(psi2[:,32]))
                global psi_vert3 = abs.(rfft(psi3[:,32]))

                # defining initial diagnostics
                E = MultiLayerQG.energies(prob)
                specE = MultiLayerQG.spectralfluxes(prob)
                fluxE = MultiLayerQG.fluxes(prob)

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

                LF1 = [fluxE[1][1]]
                LF2 = [fluxE[1][2]]
                LF3 = [fluxE[1][3]]
                VF32 = [fluxE[2][1]]
                VF52 = [fluxE[2][2]]
                TF = [fluxE[3]]

            else

                global psi1_ot = cat(psi1_ot,psi1[:,Int(round(Nx/2))], dims=2)
                global psi2_ot = cat(psi2_ot,psi2[:,Int(round(Nx/2))], dims=2)
                global psi3_ot = cat(psi3_ot,psi3[:,Int(round(Nx/2))], dims=2)

                push!(t_hovm,clock.t)

                global psi_vert1 = cat(psi_vert1,abs.(rfft(psi1[:,32])),dims=2)
                global psi_vert2 = cat(psi_vert2,abs.(rfft(psi2[:,32])),dims=2)
                global psi_vert3 = cat(psi_vert3,abs.(rfft(psi3[:,32])),dims=2)
            
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
    
                # push flux terms
                push!(LF1,fluxE[1][1])
                push!(LF2,fluxE[1][2])
                push!(LF3,fluxE[1][3])
                push!(VF32,fluxE[2][1])
                push!(VF52,fluxE[2][2])
                push!(TF,fluxE[3])
            end



            # finding vertical structure of instability
            psi_vert = [maximum(psi_vert1[:,end])
                        maximum(psi_vert2[:,end])
                        maximum(psi_vert3[:,end])]
            
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

            # save output and reset params every year
            if ((tiempo - yr_cnt*365*86400) > 0)
                # check to see if model has reached s.s.
                if ss_yr==true
                    ss_yr_cnt += 1
                end

                if KE1[end]/KE1[1] < KE_thresh && ss_yr==false
                    ss_yr==true
                    ss_yr_cnt = 1
                end


                if save_output
                    # saving output data
                    println("Saving output data to CSV for year: "*string(yr_cnt))
            
                    csv_name = "./data/threelayer_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)* ".csv"
            
                    # should I add streamfunction or PV here?? How would I use them?
                    csv_data = Dict("t" => tiempo, "CV32" => CV32, "CV52" => CV52, "KE1" => KE1, "KE2" => KE2, "KE3" => KE3, "Nz" => nlayers, "L" => L, "H" => H, "rho" => rho, "U" => U,
                                    "dt" => dt, "F_profile" => params.F, "beta" => β, "h0" => h0, "kt" => kt,
                                    "psi1_ot" => psi1_ot, "psi2_ot" => psi2_ot,
                                    "psi3_ot" => psi3_ot, "t_hovm" => t_hovm, "k_growth_lsa" => k_x[:], "sigma_ls" => sigma_LS_all,"cfl_set" => cfl_glob,
                                    "rd_LSA" => rd1, "max_evec" => max_eve1,
                                    "max_eval" => max_eva1, "PE32" => PE32, "PE52" => PE52, "CT" => CT, "NL1" => NL1, "NL2" => NL2, "NL3" => NL3, "psivert1" => psi_vert1,
                                    "psivert2" => psi_vert2, "psivert3" => psi_vert3, "alpha" => alpha, "gamma" => gamma, "LF1" => LF1, "LF2" => LF2,
                                    "LF3" => LF3, "VF32" => VF32, "VF52" => VF52, "TF" => TF)
            
                    CSV.write(csv_name, csv_data)
                end

                t_hovm = nothing
            end

        end

        global j+=1
        ###############################

    end

    if calc_growth_rate==true && perform_ls==true
        plot_growth_rate(k_x[:],sigma_LS_mid,k_emp,sigma_emp,Lx,plotpath_main)
    elseif perform_ls==true
        plot_growth_rate(k_x[:],sigma_LS_mid,[],[],Lx,plotpath_main)
    end

    # a couple of random params
    H_t = f0^2 * (Lx/kt)^2 / (g*(rho[3]-rho[2])/rho[1])

    Ri = [((0.5*(H[1]+H[2])*g*(rho[2]-rho[1])/rho[1])/((U[1]-U[2])^2))
            ((0.5*(H[2]+H[3])*g*(rho[3]-rho[2])/rho[1])/((U[2]-U[3])^2))]

    Bu_11 = (H[1]*g*(rho[2]-rho[1])/rho[1])/f0^2 * Lx^-2
    Bu_12 = (H[2]*g*(rho[2]-rho[1])/rho[1])/f0^2 * Lx^-2
    Bu_22 = (H[2]*g*(rho[3]-rho[2])/rho[1])/f0^2 * Lx^-2
    Bu_23 = (H[3]*g*(rho[3]-rho[2])/rho[1])/f0^2 * Lx^-2

    Bu = [Bu_11, Bu_12, Bu_22, Bu_23]

    inv2_Bu = sqrt.(Bu).^-1

    # save model results, if asked
    
    if save_output
        println("Saving output data to CSV")

        csv_name = "./data/threelayer_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) * "_final.csv"
        # ψ₁, ψ₂ = vars.ψ[:, :, 1], vars.ψ[:, :, 2]

        # should I add streamfunction or PV here?? How would I use them?
        csv_data = Dict("t" => tiempo, "CV32" => CV32, "CV52" => CV52, "KE1" => KE1, "KE2" => KE2, "KE3" => KE3, "Nz" => nlayers, "L" => L, "H" => H, "rho" => rho, "U" => U,
                        "dt" => dt, "F_profile" => params.F, "beta" => β, "h0" => h0, "kt" => kt,"k_growth_lsa" => k_x[:], "sigma_ls" => sigma_LS_all,"cfl_set" => cfl_glob, "H_T_scale" => H_t,
                        "Ri" => Ri, "Bu" => Bu, "inv_sqrt_Bu" => inv2_Bu, "rd_LSA" => rd1, "max_evec" => max_eve1,
                        "max_eval" => max_eva1, "PE32" => PE32, "PE52" => PE52, "CT" => CT, "NL1" => NL1, "NL2" => NL2, "NL3" => NL3, "psivert1" => psi_vert1,
                        "psivert2" => psi_vert2, "psivert3" => psi_vert3, "alpha" => alpha, "gamma" => gamma, "LF1" => LF1, "LF2" => LF2,
                        "LF3" => LF3, "VF32" => VF32, "VF52" => VF52, "TF" => TF)

        CSV.write(csv_name, csv_data)
    end
end     # end for loop thru gammas, etc.
