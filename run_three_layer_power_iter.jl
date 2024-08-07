# This is the kernel for running QG models using GeophysicalFlows.jl.
# You pass in all relevant params, and this script executes n_iter
# iterations of the power renormalization method.
# Still an open question: How do we know what initial conditions are
# appropriate, given that the condition for reiteration is 100x I.C.?
# ANSWER: Run some models to steady-state../choose I.C.s that are a bit
# less than 2 orders of magnitude than S.S.

all_files = readdir(data_dir)

if topo_type=="y_slope"
    f_name(gam,alp,h,k) = "threelayer_"*run_type*"_gamma"*string(gam)*"_alpha"*string(alp)*"_h0"* string(round(h*Lx,digits=9))*"_kt"* string(Int(k)) *"_res" * string(Int(Nx)) * ".jld2"
else
    f_name(gam,alp,h,k) = "threelayer_"*run_type*"_gamma"*string(gam)*"_alpha"*string(alp)*"_h0"* string(Int(h))*"_kt"* string(Int(k)) *"_res" * string(Int(Nx)) * ".jld2"
end

for gamma=gammas; for (i,alpha)=enumerate(alphas); for h0=h0s; for kt=kts

    fil = findall(x->x==f_name(gamma,alpha,h0,kt),all_files)
    if isempty(fil)==true

        nsubs = nsubs_all[i]

        # change variable params
        U[1] = U[2] + alpha*(H32/H52)*(U[2]-U[3])

        # rho1 = rho[2] - (abs(U[2]-U[1])*(rho[3]-rho[2]))/abs(U[2]-U[3])/gamma
        rho1 = rho[2] - gamma * (rho[3] - rho[2]) * (H32/H52)

        rho[1] = ρ[1] = rho1

        function topo_rand(h_rms, kt, Lx, Nx)
            # Borrowed from Matt Pudig..will change when I start looking at random topo.
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
                        eta_out[i,j] = 0. # (f0/H[end]) * ((h0*Lx) * ((j-Ny/2)/Ny))
                    elseif type=="x_slope"
                        eta_out[i,j] = 0. # (f0/H[end]) * ((h0*Ly) * ((i-Nx/2)/Nx))
                    end
                end
            end

            if type=="rand"
                eta_out = (f0/H[end]) * topo_rand(h0,kt,Lx,Nx)
            end

            return eta_out
        end

        topographic_pv_gradient = (0., 0.)
        aliased_fraction=1/3; T=Float64;
        grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction, T)
        if topo_type=="eggshell"
            eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"eggshell")
        elseif topo_type=="sinusoid"
            eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sinusoid")
        elseif topo_type=="y_slope"
            # eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"y_slope")
            eta = nothing
            topographic_pv_gradient = (0., h0*(f0/H[end]))
        elseif topo_type=="x_slope"
            eta = nothing
            topographic_pv_gradient = (h0*(f0/H[end]), 0.)
        elseif topo_type=="rand"
            eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"rand")
        else
            eta = 0.
        end

        println("Done with topo.")

        # Perform linear stability analysis, if asked
        if perform_ls==true
            # load module
            println("Starting LSA.")

            # perform stability analysis
            # eta=0;
            eve1,eva1,max_eve1,max_eve_phase1,max_eva1,k_x,k_y,qx1,qy1,rd1 = LinStab.lin_stab(U,V,H,beta,0.0,Nx,Ny,rho,f0,g,Float64(Lx),Float64(Ly))
            sigma_LS_all = Array(imag(eva1))
            sigma_LS_mid = sigma_LS_all[:,round(Int,Nx/2)]

            # plot_Qy(H,qy1',plotpath_main)

            global bulk_Bu = (real(rd1[2])/Lx)^2
        end

        if adjust_domain_width==true
            global L = global Lx = global Ly = dom_width_factor * real(rd1[2])

            if topo_type=="eggshell"
                eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"eggshell")
            elseif topo_type=="sinusoid"
                eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sinusoid")
            elseif topo_type=="y_slope"
                eta = nothing
                topographic_pv_gradient = (0., h0*(f0/H[end]))
            elseif topo_type=="x_slope"
                eta = nothing
                topographic_pv_gradient = (h0*(f0/H[end]), 0.)
            elseif topo_type=="rand"
                eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"rand")
            else
                eta = 0.
            end

            global dx = L/Nx
            global dt = dx*cfl_glob/U[1]/5

        end

        # define the model problem
        b = @. (g/rho0)*(rho0-rho)

        prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, H, g, ρ, U, eta=eta, topographic_pv_gradient,
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
            plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/main/"
            plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi/"
            plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi_vert/"
            plotpath_anim = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/anim/"
        elseif linear
            plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/main/"
            plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi/"
            plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
        else
            plotpath_main = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/main/"
            plotpath_psi  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi/"
            plotpath_psi_vert  = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi_vert/"
            plotpath_anim = "./figs/plots_3layer"*"_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/anim/"
        end
        plotname = "snapshots"
        filename = joinpath(filepath, "2layer.jld2")

        include("./plotting_functions.jl")

        # file management
        if isfile(filename); rm(filename); end
        if !isdir(plotpath_main); mkpath(plotpath_main); end
        if !isdir(plotpath_anim); mkpath(plotpath_anim); end
        if !isdir(plotpath_psi); mkpath(plotpath_psi); end
        if !isdir(plotpath_psi_vert); mkpath(plotpath_psi_vert); end

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

        ED = [E[3][1]]
        BD1 = [E[3][2][1]]
        BD2 = [E[3][2][2]]
        BD3 = [E[3][2][3]]

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

        EBT_pct = [0.]
        BC1_pct = [0.]
        BC2_pct = [0.]

        # initial KE of upper layer, for renormalization
        KE1_0 = E[1][1]

        println("Done with LSA.")

        if LS_criterion==true && maximum(sigma_LS_all)==0.0
            println("Skipping run for gamma = "*string(gamma)*", alpha = "*string(alpha)*", h0 = "* string(round(h0*Lx,digits=9))*", kt = "* string(Int(kt))*" because (LSA) stable.")
        else
            # trying to run the model now
            startwalltime = time()

            global cyc = 0
            global ell, j = 1, 0
            global t_hovm = nothing

            while cyc<cycles

                global ell, j 
                global psi1_ot, psi2_ot, psi3_ot

                stepforward!(prob)
                MultiLayerQG.updatevars!(prob)
                local E = MultiLayerQG.energies(prob)
                local specE = MultiLayerQG.spectralfluxes(prob)
                local fluxE = MultiLayerQG.fluxes(prob)

                if j % nsubs == 0

                    # updating variables to plot 
                    if psi_hovm  && cyc==(cycles-1)

                        psi1 = vars.ψ[:, :, 1]
                        psi2 = vars.ψ[:, :, 2]
                        psi3 = vars.ψ[:, :, 3]

                        if isnothing(t_hovm)
                            global psi1_ot = psi1[:,Int(round(Nx/2))]
                            global psi2_ot = psi2[:,Int(round(Nx/2))]
                            global psi3_ot = psi3[:,Int(round(Nx/2))]
                            global t_hovm = Array([clock.t])

                            # global psi_vert1 = abs.(rfft(psi1[:,32]))
                            # global psi_vert2 = abs.(rfft(psi2[:,32]))
                            # global psi_vert3 = abs.(rfft(psi3[:,32]))
                        else
                            global psi1_ot = cat(psi1_ot,psi1[:,Int(round(Nx/2))], dims=2)
                            global psi2_ot = cat(psi2_ot,psi2[:,Int(round(Nx/2))], dims=2)
                            global psi3_ot = cat(psi3_ot,psi3[:,Int(round(Nx/2))], dims=2)

                            push!(t_hovm,clock.t)

                            # global psi_vert1 = cat(psi_vert1,abs.(rfft(psi1[:,32])),dims=2)
                            # global psi_vert2 = cat(psi_vert2,abs.(rfft(psi2[:,32])),dims=2)
                            # global psi_vert3 = cat(psi_vert3,abs.(rfft(psi3[:,32])),dims=2)
                        end
                    end
                    
                    # push to time
                    push!(tiempo,clock.t)

                    # push to layer-wise KE
                    push!(KE1,E[1][1]/H[1])
                    push!(KE2,E[1][2]/H[2])
                    push!(KE3,E[1][3]/H[3])

                    # push to interface potential energies
                    push!(PE32,E[2][1]/((H[1]+H[2])/2))
                    push!(PE52,E[2][2]/((H[2]+H[3])/2))

                    # # push to bottom drag and biharmonic dissipation
                    # push!(ED,E[3][1])
                    # push!(BD1,E[3][2][1])
                    # push!(BD2,E[3][2][2])
                    # push!(BD3,E[3][2][3])

                    # # push to spectral flux terms
                    # push!(CV32,[specE[2][:,:,1]])
                    # push!(CV52,[specE[2][:,:,2]])
                    # push!(CL1,[specE[1][:,1]])
                    # push!(CT,[specE[3]])
                    # push!(NL1,[specE[4][:,:,1]])
                    # push!(NL2,[specE[4][:,:,2]])
                    # push!(NL3,[specE[4][:,:,3]])

                    # # push flux terms
                    # push!(LF1,fluxE[1][1])
                    # push!(LF2,fluxE[1][2])
                    # push!(LF3,fluxE[1][3])
                    # push!(VF32,fluxE[2][1])
                    # push!(VF52,fluxE[2][2])
                    # push!(TF,fluxE[3])

                    # finding vertical structure of instability
                    # psi_vert = [maximum(psi_vert1[:,end])
                    #             maximum(psi_vert2[:,end])
                    #             maximum(psi_vert3[:,end])]
                    
                    # psi_vert = psi_vert./maximum(psi_vert)
                
                    # plotting stuff
                    global plot_model
                    if plot_model==true
                        # plot_three_layer(tiempo,[KE1 KE2 KE3],[CV32[ell+1] CV52[ell+1] CL1[ell+1] CT[ell+1] NL1[ell+1] NL2[ell+1] NL3[ell+1]],vars.q,vars.v,grid,kt,h0,plotpath_main,plotname,ell) 

                        # Should I also plot the LS most unstable vert structure and the model output for most unstable wavenumber?
                        # plot_unstable_vert(H,max_eve1,psi_vert,plotpath_psi,plotname,ell)

                        # plot_layerwise_spectra(grid.kr*Lx/(2*pi),[psi_vert1[:,end] psi_vert2[:,end] psi_vert3[:,end]],plotpath_psi_vert,plotname,ell)
                    end

                    if plot_box_bool==true
                        psi1_full = vars.ψ[:, :, 1]
                        psi2_full = vars.ψ[:, :, 2]
                        psi3_full = vars.ψ[:, :, 3]
                        
                        plot_box(psi1_full,psi2_full,psi3_full,Lx,Nx,h0,plotpath_anim,plotname,ell)
                    end

                    if compute_modes_bool==true
                        modes,rd,kxr2,kxr,S = vert_modes(rho,h0,g,f0,Lx,Nx)
                    
                        psi1_full = mean(vars.ψ[:,:,1],dims=2)
                        psi2_full = mean(vars.ψ[:,:,2],dims=2)
                        psi3_full = mean(vars.ψ[:,:,3],dims=2)
                    
                        psi_slice = [psi1_full'; psi2_full'; psi3_full'];
                    
                        mode_amp_field = project_onto_modes(modes,psi_slice)
                    
                        spec = modal_nrg_spec(kxr2,S,mode_amp_field)
                    
                        modal_amp = spec_integration(kxr,spec)
                    
                        push!(EBT_pct, modal_amp[1]/sum(modal_amp))
                        push!(BC1_pct, modal_amp[2]/sum(modal_amp))
                        push!(BC2_pct, modal_amp[3]/sum(modal_amp))
                    end

                    # increase counter
                    global ell+=1
                
                    # reading out stats
                    cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
                
                    log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1: %.3e, KE2: %.3e, PE: %.3e, walltime: %.2f min",
                                    clock.step, clock.t, cfl, E[1][1], E[1][2], E[2][1], (time()-startwalltime)/60)
                    
                    # println(log)
                

                end

                global j+=1

                R = KE1_0/E[1][1] # renormalize (KE_{ref}/KE).
                if R<Rthresh

                    if cyc<cycles-1
                        global t_hovm = nothing
                    end

                    if cyc == cycles-1
                        # updating variables to plot
                        global psi1 = vars.ψ[:, :, 1]
                        global psi2 = vars.ψ[:, :, 2]
                        global psi3 = vars.ψ[:, :, 3]
                    
                        if psi_hovm  # && cyc==(cycles-1) 
                            if isnothing(t_hovm)
                                global psi1_ot = psi1[:,Int(round(Nx/2))]
                                global psi2_ot = psi2[:,Int(round(Nx/2))]
                                global psi3_ot = psi3[:,Int(round(Nx/2))]
                                global t_hovm = Array([clock.t])

                                # global psi_vert1 = abs.(rfft(psi1[:,32]))
                                # global psi_vert2 = abs.(rfft(psi2[:,32]))
                                # global psi_vert3 = abs.(rfft(psi3[:,32]))
                            else
                                global psi1_ot = cat(psi1_ot,psi1[:,Int(round(Nx/2))], dims=2)
                                global psi2_ot = cat(psi2_ot,psi2[:,Int(round(Nx/2))], dims=2)
                                global psi3_ot = cat(psi3_ot,psi3[:,Int(round(Nx/2))], dims=2)

                                push!(t_hovm,clock.t)

                                # global psi_vert1 = cat(psi_vert1,abs.(rfft(psi1[:,32])),dims=2)
                                # global psi_vert2 = cat(psi_vert2,abs.(rfft(psi2[:,32])),dims=2)
                                # global psi_vert3 = cat(psi_vert3,abs.(rfft(psi3[:,32])),dims=2)
                            end
                        end

                        # q1 = transpose(vars.q[:, :, 1])
                        # q2 = transpose(vars.q[:, :, 2])
                        # q3 = transpose(vars.q[:, :, 3])
                        
                        # push to time
                        push!(tiempo,clock.t)

                        # push to layer-wise KE
                        push!(KE1,E[1][1]/H[1])
                        push!(KE2,E[1][2]/H[2])
                        push!(KE3,E[1][3]/H[3])

                        # push to interface potential energies
                        push!(PE32,E[2][1]/((H[1]+H[2])/2))
                        push!(PE52,E[2][2]/((H[2]+H[3])/2))

                        # push to bottom drag and biharmonic dissipation
                        push!(ED,E[3][1])
                        push!(BD1,E[3][2][1])
                        push!(BD2,E[3][2][2])
                        push!(BD3,E[3][2][3])

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

                        # # finding vertical structure of instability
                        # psi_vert = [maximum(psi_vert1[:,end])
                        #             maximum(psi_vert2[:,end])
                        #             maximum(psi_vert3[:,end])]
                        
                        # psi_vert = psi_vert./maximum(psi_vert)

                    else

                        # then reset stuff
                        println("Renormalization cycle done for gamma = "*string(gamma)*", alpha = "*string(alpha)*", h0 = "* string(round(h0*Lx,digits=9))*", kt = "* string(Int(kt))*".")

                        MultiLayerQG.set_q!(prob, vars.q * sqrt(R))

                    end

                    global cyc += 1

                    println("")
                    println("Completed renormalization cycle", " ", cyc, "/", cycles)
                    println("")

                end

            end

            # Get growth rate from exponential fit to upper-layer KE time series.
            # first calculate CSP criterion
            global psi1_ot, psi2_ot, psi3_ot
            global psi_ot = [[psi1_ot[:,1:end-1]] [psi2_ot[:,1:end-1]] [psi3_ot[:,1:end-1]]]
            global cr, cr_dopp = calc_phase_speeds(psi_ot,t_hovm[1:end-1],qy1,U,Lx,Nx)
            global cr1_dopp = cr_dopp[1]; global cr2_dopp = cr_dopp[2]; global cr3_dopp = cr_dopp[3];

            psi_now = [maximum(abs.(psi1_ot[:,end])), maximum(abs.(psi2_ot[:,end])), maximum(abs.(psi3_ot[:,end]))]
            global csp_crit, csp_terms = calc_csp_crit(cr,qy1,U,psi_now,length(H))
            if calc_growth_rate==true
                sigma_emp = LinStab.calc_growth(tiempo[1:end-1], [KE1[1:end-1] KE2[1:end-1] KE3[1:end-1] PE32[1:end-1] PE52[1:end-1]])
                sigma_emp_KE1, sigma_emp_KE2, sigma_emp_KE3 = sigma_emp[1], sigma_emp[2], sigma_emp[3]
                sigma_emp_PE32, sigma_emp_PE52 = sigma_emp[4], sigma_emp[5]
                a = findmax(CV32[end][1])
                k_emp = grid.kr[a[2][1]]
            end

            if calc_growth_rate==true && perform_ls==true
                # plot_growth_rate(k_x[:],sigma_LS_mid,k_emp,sigma_emp_KE1,Lx,plotpath_main)
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

                global data_dir

                if topo_type=="y_slope"
                    csv_name = data_dir*"/threelayer_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) * ".jld2"
                else
                    csv_name = data_dir*"/threelayer_"*run_type*"_gamma"*string(gamma)*"_alpha"*string(alpha)*"_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) * ".jld2"
                end
                
                # should I add streamfunction or PV here?? How would I use them?
                csv_data = Dict("t" => tiempo, "CV32" => CV32, "CV52" => CV52, "KE1" => KE1, "KE2" => KE2, "KE3" => KE3, "Nz" => nlayers, "L" => L, "H" => H, "rho" => rho, "U" => U,
                                "dt" => dt, "F_profile" => params.F, "beta" => β, "h0" => h0, "kt" => kt, "sigma_emp_KE1" => sigma_emp_KE1, "sigma_emp_KE2" => sigma_emp_KE2,
                                "sigma_emp_KE3" => sigma_emp_KE3, "sigma_emp_PE32" => sigma_emp_PE32, "sigma_emp_PE52" => sigma_emp_PE52, "psi1_ot" => psi1_ot, "psi2_ot" => psi2_ot,
                                "psi3_ot" => psi3_ot, "t_hovm" => t_hovm, "k_growth_emp" => k_emp, "k_growth_lsa" => k_x[:], "sigma_ls" => sigma_LS_all,"cfl_set" => cfl_glob, "H_T_scale" => H_t,
                                "Ri" => Ri, "Bu" => Bu, "inv_sqrt_Bu" => inv2_Bu, "csp_crit" => csp_crit, "csp_terms" => csp_terms, "rd_LSA" => rd1, "max_evec" => max_eve1,
                                "max_eval" => max_eva1, "PE32" => PE32, "PE52" => PE52, "CT" => CT, "NL1" => NL1, "NL2" => NL2, "NL3" => NL3, "Qy" => params.Qy,
                                "alpha" => alpha, "gamma" => gamma, "cr" => cr, "cr_Dopp" => cr_dopp, "LF1" => LF1, "LF2" => LF2, "EBT_pct" => EBT_pct, "BC1_pct" => BC1_pct, "BC2_pct" => BC2_pct,
                                "LF3" => LF3, "VF32" => VF32, "VF52" => VF52, "TF" => TF, "Ekman_drag" => ED, "biharmonic_diss_1" => BD1, "biharmonic_diss_2" => BD2,
                                "biharmonic_diss_3" => BD3, "eta" => eta, "psi1_full" => psi1, "psi2_full" => psi2, "psi3_full" => psi3, "nsubs" => nsubs)

                jldsave(csv_name; csv_data)
            end
        end
    else
        println("File for gamma: "*string(gamma)*", alpha: "*string(alpha)*", h0: "*string(h0)*", kt: "*string(kt)*" already exists.")
    end

end; end; end; end

