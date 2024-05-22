# This is the kernel for running QG models using GeophysicalFlows.jl.

for U1=U1s; for rho1=rho1s

    # change variable params
    U[1] = U1

    rho[1] = ρ[1] = rho1

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
    else
        eta = 0.
    end

    b = zeros(nlayers)
    b[1] = (g/rho0)*(rho0-rho[1])
    b[2] = (g/rho0)*(rho0-rho[2])

    # define the model problem
    prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, H, b, U, nν, ν,
    μ, β, dt, stepper, linear, aliased_fraction=1/3)

    sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid
    x, y = grid.x, grid.y

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
        plotpath_main = "./figs/plots_2layer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_linear_res" * string(Int(Nx)) *"/main/"
    else
        plotpath_main = "./figs/plots_2layer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_res" * string(Int(Nx)) *"/main/"
    end

    plotname = "snapshots"

    include("./plotting_functions.jl")

    # file management
    if !isdir(plotpath_main); mkpath(plotpath_main); end

    # trying to run the model now
    startwalltime = time()

    global ell, j = 1, 0
    global t_yrly = nothing
    global yr_cnt = 1
    global ss_yr = false
    global ss_yr_cnt = 0

    while ss_yr_cnt < ss_yr_max
        global ell, j 

        ##########################
        stepforward!(prob)
        MultiLayerQG.updatevars!(prob)
        local E = MultiLayerQG.energies(prob)

        if j % nsubs == 0

            if isnothing(t_yrly)
                global psi1_ot = vars.ψ[:, :, 1]
                global psi2_ot = vars.ψ[:, :, 2]

                global q1_ot = vars.q[:,:,1]
                global q2_ot = vars.q[:,:,2]

                global t_yrly = Array([clock.t])

                # defining initial diagnostics
                E = MultiLayerQG.energies(prob)

                # variables for plotting, that will be pushed to
                global KE1 = [E[1][1]]/H[1]
                global KE2 = [E[1][2]]/H[2]

                global PE32 = [E[2][1]]/((H[1]+H[2])/2)

            else

                global psi1_ot = cat(psi1_ot, vars.ψ[:, :, 1], dims=3)
                global psi2_ot = cat(psi2_ot, vars.ψ[:, :, 2], dims=3)

                global q1_ot = cat(q1_ot, vars.q[:, :, 1], dims=3)
                global q2_ot = cat(q2_ot, vars.q[:, :, 2], dims=3)
                
                # defining initial diagnostics
                E = MultiLayerQG.energies(prob)

                # push to time
                push!(t_yrly,clock.t)

                # push to layer-wise KE
                push!(KE1,E[1][1]/H[1])
                push!(KE2,E[1][2]/H[2])
    
                # push to interface potential energies
                push!(PE32,E[2][1]/((H[1]+H[2])/2))

                # plotting stuff
                global plot_model
                if plot_model==true

                end

                # increase counter
                global ell+=1

            end
        
            # reading out stats
            cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
        
            log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE1: %.3e, KE2: %.3e, PE: %.3e, walltime: %.2f min",
                            clock.step, clock.t, cfl, E[1][1], E[1][2], E[2][1], (time()-startwalltime)/60)
            
            println(log)

            # save output and reset params every year
            if ((t_yrly[end] - yr_cnt*365*86400) > 0)
                yr_cnt += 1
                ell = 0
                # check to see if model has reached s.s.
                if ss_yr==true
                    ss_yr_cnt += 1
                end

                if KE1[end]/KE1[1] < KE_thresh && ss_yr==false
                    ss_yr = true
                    ss_yr_cnt = 1
                end

                if isnan(KE1[end])
                    ss_yr_cnt = ss_yr_max
                end

                if save_output
                    # saving output data            
		            println("Saving annual data for year: "*string(yr_cnt))
                    if linear
                        jld_name = data_dir*"/twolayer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_linear_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld2"
                    else
                        jld_name = data_dir*"/twolayer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_NL_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld2"
                    end

		        println("Saving output data to JLD to: "*jld_name)

                jld_data = Dict("t" => t_yrly, "KE1" => KE1, "KE2" => KE2, "Nz" => nlayers,
                                "L" => L, "H" => H, "rho" => rho, "U" => U,
                                "dt" => dt, "beta" => β,
                                "psi1_ot" => Array(psi1_ot), "psi2_ot" => Array(psi2_ot),
                                "q1_ot" => Array(q1_ot), "q2_ot" => Array(q2_ot),
                                "cfl_set" => cfl_glob, "PE32" => PE32)
            
                jldsave(jld_name; jld_data)

                end

                t_yrly = nothing
            end

        end

        global j+=1
        ###############################

    end

    # save model results, if asked
    
    if save_output
        println("Saving output data to JLD")
        if linear
            jld_name = data_dir*"/twolayer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_linear_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld2"
        else
            jld_name = data_dir*"/twolayer_U_"*string(round(U[1],sigdigits=5))*"_rho_"*string(round(rho[1],sigdigits=5))*"_NL_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld2"
        end
        
        # should I add streamfunction or PV here?? How would I use them?
        jld_data = Dict("t" => t_yrly, "KE1" => KE1, "KE2" => KE2,
                        "Nz" => nlayers, "L" => L, "H" => H, "rho1s" => rho1s, "U1s" => U1s,
                        "dt" => dt, "beta" => β, "cfl_set" => cfl_glob, "PE32" => PE32)

        jldsave(jld_name; jld_data)
    end
end; end    # end for loop thru U1s, rho1s
