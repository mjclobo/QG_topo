# This is a module of functions used for running QG models using GeophysicalFlows.jl.
# 

####################################################################################
## User input parameters
####################################################################################

@with_kw struct mod_params
    data_dir::String = "./"
    Nz::Int64 = 2
    Nx::Int64 = 562
    Ny::Int64 = 562
    Lx::Int64 = 1000.e3
    Ly::Int64 = 1000.e3
    Ld::Float64 = 25.e3
    f0::Float64 = 8.e-5
    H::Array{Float64} = ones(Nz) * (4000/Nz)
    g::Float64 = 9.81
    rho0::Float64 = 1025.
    rho::Array{Float64} = zeros(Nz)
    strat_str::String = "uni_strat"
    b::Array{Float64} = (g/rho0) .* (rho0 .- rho)
    shear_str::String = "uni_shear"
    U::Array{Float64} = ones(Nz)
    μ::Float64 = 0.
    κ::Float64 = 0.
    nν::Int64 = 4
    ν::Float64 = 0.
    q0_mag::Float64 = 1.e-7
    eta::Any = zeros(Nx,Ny) # eta::Array{Float64} = zeros(Nx,Ny)
    topographic_pv_gradient::Tuple{Float64, Float64} = (0., 0.)
    topo_type::String = "flat"
    h0::Tuple{Float64, Float64} = (0., 0.)
    kt::Float64 = 0.
    β::Float64 = 0.
    dt::Float64 = 1200.
    stepper::String = "ETDRK4"
    linear::Bool = false
    dev::Any = CPU()
    restart_bool::Bool = false
    restart_yr::Int = 0
    yr_increment::Float64 = 1.0
    ss_yr_max::Int = 100
    nsubs::Int64 = round(Int64, 5*(Ld / (U[1]/2))/dt)  # save psi field every 5 eddy periods
end

####################################################################################
## File naming convention
####################################################################################

function jld_name(model_params, yr_cnt)

    @unpack_mod_params model_params

    # topo string
    if topo_type=="y_slope"
        h_str = "_" * topo_type * "_h0_" * (@sprintf "%.3E" h0[1])
    elseif topo_type=="rough"
        h_str = "_" * topo_type * "_hrms_" * string(round(h0[2]))*"_kt" * string(round(kt)) 
    elseif topo_type=="rand_slope"
        h_str = "_" * topo_type * "_hrms_"* string(round(h0[2]))*"_kt" * string(round(kt)) 
    elseif topo_type=="sin_sin"
        h_str = "_" * topo_type * "_h0_"* string(Int(h0[1]))  
    elseif topo_type=="flat"
        h_str = "_" * topo_type
    end

    # drag strings
    kap_str = @sprintf "%.3E" κ
    mu_str = @sprintf "%.3E" μ

    nu_str = @sprintf "%.3E" ν
    hv_str = "_nu" * nu_str

    drag_str = "_mu" * mu_str * "_kappa" * kap_str

    # misc.
    beta_str = "_beta" * (@sprintf "%.3E" β) * "_"

    res_str =  "_res" * string(Int(Nx))

    yr_str = "_yr" * string(yr_cnt)

    # geometry
    L_str = "_L_2pi" * (@sprintf "%.3E" Lx)

    thick_str = "layer" * string(Nz) * "_Htot_" * string(sum(H))
    
    return "/" * thick_str * L_str * h_str * beta_str * shear_str * "_" * strat_str * drag_str * hv_str * res_str * yr_str * ".jld"
end

####################################################################################
## Topographic stuff
####################################################################################
# must set topographic_pv_gradient (tuple of the large-scale slope)
# and eta (periodic topography)

function topo_rough(h_rms, kt, Lx, Nx)
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
                eta_out[i,j] = (f0/H[end]) * ((h0*Lx) * ((j-Ny/2)/Ny))
            elseif type=="sin_sin"
                eta_out[i,j] = (f0/H[end]) * (h0[1] * sin(2*pi*kt[1]*y[j]/Ly) + h0[2] * sin(2*pi*kt[2]*y[j]/Ly))
            end
        end
    end

    if type=="rand"
        eta_out = (f0/H[end]) * topo_rough(h0,kt,Lx,Nx)
    end

    return eta_out
end

function define_topo(model_params)

    @unpack_mod_params model_params

    topographic_pv_gradient = (0., 0.)

    T=Float64;
    grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction=1/3, T)

    if topo_type=="eggshell"
        eta = topographicPV(grid_topo,h0[2],kt,Lx,Ly,f0,H,"eggshell")
    elseif topo_type=="sinusoid"
        eta = topographicPV(grid_topo,h0[2],kt,Lx,Ly,f0,H,"sinusoid")
    elseif topo_type=="sin_sin"
        eta = topographicPV(grid_topo,h0[2],kt,Lx,Ly,f0,H,"sin_sin")
    elseif topo_type=="y_slope"
        # eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"y_slope")
        eta = nothing
        topographic_pv_gradient = (0., h0[1]*(f0/H[end]))
    elseif topo_type=="rand_slope"
        eta = topographicPV(grid_topo,h0[2],kt,Lx,Ly,f0,H,"rand")
        topographic_pv_gradient = (0., h0[1]*(f0/H[end]))
    elseif topo_type=="rand_flat"
        eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"rand")
        topographic_pv_gradient = (0., 0.)
    else
        eta = nothing
    end

    return topographic_pv_gradient, eta
end

####################################################################################
## Functions to initialize model
####################################################################################

function set_initial_conditions(prob, prob_filt, model_params)
    
    @unpack_mod_params model_params
    
    # setting initial conditions
    grid = prob.grid;
    
    if restart_bool==false
        seed!(1234) # reset of the random number generator for reproducibility
        q₀  = q0_mag * device_array(dev)(randn((Nx, Ny, Nz)))
        q₀h = prob_filt.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
        q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

        MultiLayerQG.set_q!(prob, q₀)

    elseif restart_bool==true

        r_file_name = jld_name(model_params, restart_yr)

        a = load(data_dir * r_file_name)

        ψ = a["jld_data"]["psi_ot"][:,:,:,1]

        MultiLayerQG.set_ψ!(prob.sol, prob.params, prob.vars, prob.grid, ψ)   # this also sets q!!!

    end

    return prob
end


function initialize_model(model_params)

    @unpack_mod_params model_params

    nlayers = Nz

    # this just defines wavenumber filter, in case model doesn't use one
    stepper2 = "FilteredRK4"
    prob_filt = MultiLayerQG.Problem(nlayers, dev; nx=Nx, Lx=Lx, f₀=f0, H, b, U, nν, ν, eta, topographic_pv_gradient,
        μ, κ, β, dt, stepper=stepper2, linear, aliased_fraction=1/3)

    prob = MultiLayerQG.Problem(nlayers, dev; nx=Nx, Lx=Lx, f₀=f0, H, b, U, nν, ν, eta, topographic_pv_gradient,
        μ, κ, β, dt, stepper, linear, aliased_fraction=1/3)

    # sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid
    # x, y = grid.x, grid.y

    return prob, prob_filt
end

####################################################################################
## Run model
####################################################################################
# trying to run the model now

function run_model(prob, model_params)

    @unpack_mod_params model_params

    sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid

    startwalltime = time()

    global j = 0
    global t_yrly = nothing
    global yr_cnt = restart_yr

    while yr_cnt < ss_yr_max
        global j
        global prob
        
        stepforward!(prob)
        MultiLayerQG.updatevars!(prob)

        if j % nsubs == 0

            if isnothing(t_yrly)
                global psi_ot = deepcopy(vars.ψ);

                global t_yrly = Array([clock.t])
            else
                global psi_ot = cat(psi_ot, vars.ψ, dims=4)
                
                # push to time
                push!(t_yrly,clock.t)
            end
        
            # reading out stats
            cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
        
            log = @sprintf("step: %04d, t: %.1f day, cfl: %.2f, walltime: %.2f min, nu_star: %.2E",
                            clock.step, clock.t/3600/24, cfl, (time()-startwalltime)/60, prob.params.ν)
            
            println(log)

            # save output and reset params every year
            if ((t_yrly[end] - yr_cnt*365*86400) > 0)

                if yr_cnt==0.0

                    jld_data = Dict("t" => t_yrly, "Nz" => Nz,
                        "L" => L, "H" => H, "rho" => rho, "U" => U[:,:,1],
                        "dt" => dt, "beta" => β, "mu" => μ, "kappa" => κ,
                        "psi_ot" => Array(psi_ot),
                        "nu" => ν, "n_nu" => nν, "Ld" => Ld,
                        "Qy" => Array(params.Qy[1,1,:]), "eta" => eta,
                        "topographic_pv_gradient" => topographic_pv_gradient,
                        "stepper" => stepper)

                else

                    jld_data = Dict("t" => t_yrly,
                        "psi_ot" => Array(psi_ot))

                end

                # saving output

                save_output(vars, jld_data, model_params, yr_cnt)

                GC.gc()

            end

        end

        global j+=1
        ###############################

    end

end


####################################################################################
## Save output
####################################################################################

function save_output(vars, jld_data, model_params, yr_cnt)

    @unpack_mod_params model_params

    # saving yearly output
    println("Saving annual data for year: " * string(yr_cnt))

    if dev==GPU()
        psi1 = CUDA.@allowscalar vars.ψ[1,1,1]
    else
        psi1 = vars.ψ[1,1,1]
    end

    if isnan(psi1)
        global yr_cnt = ss_yr_max
    else
        global yr_cnt += yr_increment
    end

    file_name = jld_name(model_params, yr_cnt - yr_increment)
    
    println("Saving output data to JLD to: " * file_name)

    jldsave(data_dir * file_name; jld_data)

    global t_yrly = nothing

end


####################################################################################
## 
####################################################################################





####################################################################################
## density profile; currently doing a naive approach, should switch to
# a minimum rmse approach, though
####################################################################################


function vert_disc(strat_type, rho_top, rho_bottom, H, scale_depth)
    
    Nz = length(H)

    rho_out = zeros(Nz)

    z_prof = zeros(Nz)

    for i in range(1,Nz)
        z_prof[i] = - (sum(H[1:i-1]) + H[i]/2)
    end

    for (i,z) in enumerate(z_prof)
        rho_out[i] = vert_profile(strat_type, rho_top, rho_bottom, H, scale_depth, z)
    end

    return rho_out
end


function vert_profile(strat_type, rho_top, rho_bottom, H, scale_depth, z)
    if strat_type[1:2] == "SI"
        return rho_top + (rho_bottom - rho_top) * (1 - exp(z / scale_depth))
    elseif strat_type[1:3] == "uni"
        return rho_top - (rho_bottom - rho_top) * (z/sum(H))
    end
end


# vert_disc("SI_strat", 1026, 1027.5, 1000.0 * ones(4), 750)

####################################################################################
## Auxiliary functions
####################################################################################

function gp(rho,rho0,g)
    # g_prime = g*(rho[2]-rho[1])/rho0
    g_prime = g*(rho[2]-rho[1])/rho0
    return g_prime
end

function calc_stretching_mat(model_params; rigid_lid=true)

    @unpack_mod_params model_params

    S = zeros((Nz,Nz,))

    if rigid_lid
        alpha = 0
    else
        alpha = -f0^2/g/H[1]
    end

    S[1,1] = -f0^2/H[1]/gp(rho[1:2],rho0,g) + alpha
    S[1,2] = f0^2/H[1]/gp(rho[1:2],rho0,g)
    
    for i = 2:Nz-1
        S[i,i-1] = f0^2/H[i]/gp(rho[i-1:i],rho0,g)
        S[i,i]   = -(f0^2/H[i]/gp(rho[i-1:i],rho0,g) + f0^2/H[i]/gp(rho[i:i+1],rho0,g))
        S[i,i+1] = f0^2/H[i]/gp(rho[i:i+1],rho0,g)
    end

    S[Nz,Nz-1] = f0^2/H[Nz]/gp(rho[Nz-1:Nz],rho0,g)
    S[Nz,Nz]   = -f0^2/H[Nz]/gp(rho[Nz-1:Nz],rho0,g)

    return S
end


function calc_Ld(model_params)

    S = calc_stretching_mat(model_params)

    return abs.(eigvals(S)).^-0.5
end
