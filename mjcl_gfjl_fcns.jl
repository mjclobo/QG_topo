# This is a module of functions used for running QG models using GeophysicalFlows.jl.
# 
import Base: view
####################################################################################
## User input parameters
####################################################################################

@with_kw struct mod_params
    data_dir::String = "./"
    Nz::Int64 = 2
    Nx::Int64 = 512
    Ny::Int64 = 512
    Lx::Float64 = 1000.0e3
    Ly::Float64 = 1000.0e3
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
    dyn_nu::Bool = false
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
    restart_yr::Float64 = 0.
    pre_buoy_restart_file::Bool = false
    data_dir_pre_buoy::String = data_dir
    yr_increment::Float64 = 1.0
    ss_yr_max::Int = 100
    nsubs::Int64 = round(Int64, 5*(Ld / (U[1]/2))/dt)  # save psi field every 5 eddy periods
end

@with_kw struct diag_bools
    psi_out_bool::Bool = false
    xspace_layered_nrg::Bool = false
    kspace_layered_nrg::Bool = false
    xspace_modal_nrg::Bool = false
    kspace_modal_nrg::Bool = false
    xspace_layered_nrg_budget::Bool = false
    kspace_layered_nrg_budget::Bool = false
    xspace_modal_nrg_budget::Bool = false
    kspace_modal_nrg_budget::Bool = false
    two_layer_kspace_modal_nrg_budget_bool::Bool = false
    psi_out_bool_yrs_end::Bool = false
    only_save_last::Bool = false
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


function jld_name_pre_buoy(model_params,yr_cnt)

    @unpack_mod_params model_params
    
    drag_str="_lin_drag_"
    mu_str = @sprintf "%.2E" μ / (((U[1]) / 2) / Ld)
    
    nu_str = @sprintf "%.2E" ν / ((U[1]/2) * (Lx/2/pi)^7)
       
    return data_dir_pre_buoy*"/twolayer_L2pi_" * string(round(L/2/pi/Ld)) * "_h0"* string(round(h0[1]/S32,digits=3))* "_beta" * string(round(β * 2 * Ld^2 / U[1],digits=3)) * "_U" * string(round(U[1],digits=4)) * "_rho"* string(round(ρ[1],digits=6)) * drag_str * "mu" * mu_str * "_nu" * nu_str * "_Hr" * string(round(H[1]/H[2],sigdigits=1)) * "_res" * string(Int(Nx)) * "_yr" * string(Int(yr_cnt)) *  ".jld"
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

        if pre_buoy_restart_file==true

            r_file_name = jld_name_pre_buoy(model_params, restart_yr)

            a = load(r_file_name)

            ψ = a["jld_data"]["psi_ot"][:,:,:,1]

            MultiLayerQG.set_ψ!(prob.sol, prob.params, prob.vars, prob.grid, ψ)   # this also sets q!!!

            println("Restarting model from: " * r_file_name)

        else

            r_file_name = jld_name(model_params, restart_yr)

            a = load(data_dir * r_file_name)

            ψ = a["jld_data"]["psi_yrs_end"]

            MultiLayerQG.set_ψ!(prob.sol, prob.params, prob.vars, prob.grid, ψ)   # this also sets q!!!

            println("Restarting model from: " * data_dir * r_file_name)

        end

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

    @unpack_diag_bools diags

    @unpack_mod_params model_params

    sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)

    startwalltime = time()

    global j = 0
    prob.clock.t = restart_yr * 365.25 * 24 * 3600.
    global t_yrly = Array([prob.clock.t])
    global yr_cnt = restart_yr
    global budget_counter = 0
    global nsaves = 0
    global psi_ot = nothing

    global ph_iso    = 0.
    global ph_slices = 0.

    nterms_two_layer_modal_kspace = 14  # including residual
    nterms_two_layer_modal_xspace = 10  # only one nonlinear term (instead of 5)
    global two_layer_kspace_modal_nrgs = zeros(dev, T, (grid.nkr, nterms_two_layer_modal_kspace))
    global two_layer_xspace_modal_nrgs = zeros(dev, T, (nterms_two_layer_modal_xspace, 1))
    global two_layer_modal_length_scales = zeros(dev, T, (2, 1))
    global two_layer_vBT_scale = 0.

    len_nrg = ceil(Int, (ss_yr_max - yr_cnt + 1) * 365.25 * 24 * 3600 / prob.clock.dt / nsubs)
    global nrg_ot = zeros(dev, T, (3, len_nrg))

    while yr_cnt < ss_yr_max
        global j
        
        # global prob

        if dyn_nu==true
            rmsζ = sqrt(mean((irfft(-grid.Krsq .* prob.vars.ψh[:,:,1], grid.ny)).^2))
            global prob = @set prob.params.ν = rmsζ * dx^8
        end
        
        stepforward!(prob)
        MultiLayerQG.updatevars!(prob)

        if xspace_modal_nrg_budget==true
            # X,x = layered_xspace_budget(vars.ψ, model_params)
        end

        if kspace_modal_nrg_budget==true
            # X,x = layered_kspace_budget(vars.ψ, model_params)
        end

        if j % nsubs == 0
            global nsaves+=1

            if psi_out_bool==true
                if isnothing(psi_ot)
                    global psi_ot = deepcopy(vars.ψ);
                else
                    global psi_ot = cat(psi_ot, vars.ψ, dims=4)
                end
            end

            push!(t_yrly,prob.clock.t)

            if two_layer_kspace_modal_nrg_budget_bool==true
                
                global uBT_rms, two_layer_kspace_modal_nrgs, two_layer_xspace_modal_nrgs, two_layer_modal_length_scales = update_two_layer_kspace_modal_nrgs(prob, vars.ψ, model_params, two_layer_kspace_modal_nrgs, two_layer_xspace_modal_nrgs, two_layer_modal_length_scales)
                # global nrg_ot_here, two_layer_xspace_layer_nrgs = update_two_layered_nrg(prob, vars.ψ, model_params, two_layer_xspace_layer_nrgs) 
                global @views nrg_ot[:,nsaves] = update_two_layered_nrg(prob, vars.ψ, model_params) 
                global two_layer_vBT_scale += uBT_rms

                global ph_iso    += calc_iso_phase_shift_12(model_params, vars.ψ, prob.grid, prob.vars)
                global ph_slices += calc_zonal_phase_shift_12(model_params, vars.ψ, prob.grid, prob.vars)

                global budget_counter +=1

            end
        
            # reading out stats
            cfl = prob.clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
        
            log = @sprintf("step: %04d, t: %.1f day, cfl: %.2f, walltime: %.2f min, nu_star: %.2E",
                            clock.step, clock.t/3600/24, cfl, (time()-startwalltime)/60, prob.params.ν)
            
            println(log)

            # save output and reset params every year
            if ((t_yrly[end] - yr_cnt*365*86400) > 0)

                if yr_cnt <= restart_yr + yr_increment

                    jld_data = Dict("t" => t_yrly, "Nz" => Nz,
                        "L" => L, "H" => H, "rho" => rho, "U" => U[:,:,1],
                        "dt" => dt, "beta" => β, "mu" => μ, "kappa" => κ,
                        "nu" => ν, "n_nu" => nν, "Ld" => Ld,
                        "Qy" => Array(params.Qy[1,1,:]), "eta" => eta,
                        "topographic_pv_gradient" => topographic_pv_gradient,
                        "stepper" => stepper)

                else

                    @unpack_diag_bools diags

                    if psi_out_bool==true && two_layer_kspace_modal_nrg_budget_bool==false
                        jld_data = Dict("t" => t_yrly,
                            "psi_ot" => Array(psi_ot))
                    elseif psi_out_bool==true && two_layer_kspace_modal_nrg_budget_bool==true
                        jld_data = Dict("t" => t_yrly,
                            "psi_ot" => Array(psi_ot), "two_layer_kspace_modal_nrg_budget" => Array(two_layer_kspace_modal_nrgs ./ budget_counter),
                            "two_layer_xspace_layer_nrgs" => Array(two_layer_xspace_layer_nrgs ./ budget_counter),
                            "two_layer_vBT_scale" => Float64(two_layer_vBT_scale ./ budget_counter))
                    elseif psi_out_bool_yrs_end==true && two_layer_kspace_modal_nrg_budget_bool==true
                        jld_data = Dict("t" => t_yrly,
                            "psi_yrs_end" => Array(vars.ψ), "two_layer_kspace_modal_nrg_budget" => Array(two_layer_kspace_modal_nrgs ./ budget_counter),
                            "two_layer_xspace_modal_nrg_budget" => Array(two_layer_xspace_modal_nrgs ./ budget_counter),
                            "two_layer_xspace_nrgs_ot" => Array(nrg_ot),
                            "two_layer_vBT_scale" => Float64(two_layer_vBT_scale ./ budget_counter),
                            "ph_iso" => Float64(ph_iso / budget_counter), "ph_slices" => Float64(ph_slices / budget_counter),
                            "two_layer_modal_length_scales" => Array(two_layer_modal_length_scales ./ budget_counter))
                    elseif psi_out_bool==false && two_layer_kspace_modal_nrg_budget_bool==true
                        jld_data = Dict("t" => t_yrly,
                            "two_layer_kspace_modal_nrg_budget" => Array(two_layer_kspace_modal_nrgs ./ budget_counter),
                            "two_layer_xspace_modal_nrg_budget" => Array(two_layer_xspace_modal_nrgs ./ budget_counter),
                            "two_layer_xspace_nrgs_ot" => Array(nrg_ot),
                            "two_layer_vBT_scale" => Float64(two_layer_vBT_scale ./ budget_counter),
                            "ph_iso" => Float64(ph_iso / budget_counter), "ph_slices" => Float64(ph_slices / budget_counter),
                            "two_layer_modal_length_scales" => Array(two_layer_modal_length_scales ./ budget_counter))
                    end

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

    @unpack_diag_bools diags

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
        global yr_cnt = round(yr_cnt + yr_increment, digits=3)
    end

    file_name = jld_name(model_params, round(yr_cnt - yr_increment, digits=3))
    
    println("Saving output data to JLD to: " * file_name)

    if restart_bool==true
        if yr_cnt > round(restart_yr + yr_increment, digits=3)
            jldsave(data_dir * file_name; jld_data)
        end
    else
        jldsave(data_dir * file_name; jld_data)
    end

    if yr_cnt - 3 * yr_increment > restart_yr
        if only_save_last==true
            rm(data_dir * jld_name(model_params, round(yr_cnt - 2 * yr_increment, digits=3)))
        end
    end

    global psi_ot = nothing

end


# function set_output_dict(diags)

# end

####################################################################################
## DIAGS: x-space modal budget
####################################################################################

using Statistics

function modal_xspace_budget(model_params, vars, params, grid, sol, ψ)

    @unpack_mod_params model_params

    nlayers=Nz
    δ = params.H[1]/params.H[2]
    g′ = g * (rho[2] -   rho[1]) / rho0 # reduced gravity at the interface
    μ = params.μ
    Ld1 = (params.f₀^2 / (g′ * params.H[1]))^-0.5
    ∂yηb = params.topographic_pv_gradient[2] / (params.f₀ / params.H[end])

    vars.ψh[:,:,1] = rfft(ψ[:,:,1])
    vars.ψh[:,:,2] = rfft(ψ[:,:,2])

    ψBC = 0.5 * (ψ[:,:,1] .- ψ[:,:,2])
    ψBT = 0.5 * (ψ[:,:,1] .+ ψ[:,:,2])

    ψBCh = rfft(ψBC)
    ψBTh = rfft(ψBT)

    ∂xψBC = irfft(im * grid.kr .* ψBCh, grid.ny)
    ∂xψBT = irfft(im * grid.kr .* ψBTh, grid.ny)

    ∂yψBC = irfft(im * grid.l .* ψBCh, grid.ny)
    ∂yψBT = irfft(im * grid.l .* ψBTh, grid.ny)

    ζ = zeros(grid.nx, grid.ny, 2)   # irfft(- grid.Krsq .* ψh, grid.nx)
    ζ[:,:,1] = irfft(- grid.Krsq .* vars.ψh[:,:,1], grid.nx)
    ζ[:,:,2] = irfft(- grid.Krsq .* vars.ψh[:,:,2], grid.nx)

    ∂xζ1 = irfft(im * grid.kr .* rfft(ζ[:,:,1]), grid.ny)

    ζBC = irfft(- grid.Krsq .* ψBCh, grid.nx)  # BC
    ζBT = irfft(- grid.Krsq .* ψBTh, grid.nx)   # BT

    ∂xζBC = irfft(im * grid.kr .* rfft(ζBC), grid.ny)
    ∂yζBC = irfft(im * grid.l  .* rfft(ζBC), grid.ny)

    ∂xζBT = irfft(im * grid.kr .* rfft(ζBT), grid.ny)
    ∂yζBT = irfft(im * grid.l  .* rfft(ζBT), grid.ny)


    #     vertfluxes, dragfluxes = zeros(grid.nkr,nlayers), zeros(grid.nkr,nlayers+1)

    U₁, U₂ = view(params.U, :, :, 1), view(params.U, :, :, 2)

    S32 = CUDA.@allowscalar params.f₀ * U₁[1,1] ./ g′


    # calcing x_space_stuff
    # Linear modal transfer terms
    LT_x = @. 0.5 * ψBT * U₁ * ∂xζ1   # BC to BT

    # Topo transfer term
    TT_x = - 0.5 * ∂yηb * (params.f₀ / params.H[2]) .* ψBT .* ∂xψBC
    # k-space version: -0.5 * (params.f₀ / params.H[2]) * ∂yηb * (conj(ψBTh) * ∂xψBCh + ψBTh * conj(∂xψBCh))

    # Baroclinic production term
    BC_x = (params.f₀ / params.H[1]) * S32 .* ψBC .* ∂xψBT

    # Nonlinear transfer term
    NL_x = ψBT .* (∂xψBC .* ∂yζBC .- ∂yψBC .* ∂xζBC)  # NL BT <--> BC (EKE)

    # Drag terms
    DBT_x = @. - (μ/2) * ψBT .* (ζBC .- ζBT)     # BT drag term
    DBC_x = @. (μ/2) * ψBC .* (ζBC .- ζBT)      # BC drag term    

    # integrating in y
    LT_x = sum(LT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    TT_x = sum(TT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    BC_x = sum(BC_X) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    NL_x = sum(NL_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    DBT_x = sum(DBT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    DBC_x = sum(DBC_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    resid_x = BC_x + DBT_x + DBC_x

    # energies
    BTKE_x = 0.5 * sum(∂xψBT.^2 .+ ∂yψBT.^2) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    BCKE_x = 0.5 * sum(∂xψBC.^2 .+ ∂yψBC.^2) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    BCEAPE_x = 2 * Ld1^-2 * sum(ψBC .* ψBC)  * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    return 
end

modal_xspace_budget(model_params, prob, ψ) = modal_xspace_budget(model_params, prob.vars, prob.params, prob.grid, prob.sol, ψ)





####################################################################################
## DIAGS: modal k-space budget
####################################################################################

function update_two_layer_kspace_modal_nrgs(vars, params, grid, sol, ψ, model_params, nrgs_in, nrgs_in_x, lengths_in)
    # energies are: BTEKE, BCEKE, EAPE; CBC, DBC, DBT; Tflat, Ttopo; NLBCEAPE, NLBCEKE, NLBC2BT; NLBTEKE, NLBT2BC; resid
    # here we do not define average, just add up the budget...averaging comes later


    @unpack_mod_params model_params

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)

    rfftplan = plan_flows_rfft(A{T, 3}(undef, grid.nx, grid.ny, 1), [1, 2]; flags=FFTW.MEASURE)

    # parameters
    nlayers = 2
    δ = [params.H[1]/sum(params.H), params.H[2]/sum(params.H)]
    g′ = g * (rho[2] -   rho[1]) / rho0 # reduced gravity at each interface
    Ld1 = (f0^2 / (g′ * params.H[1]))^-0.5
    ∂yηb = params.topographic_pv_gradient[2] / (f0 / params.H[end])
    
    # assigning basic variables
    ψBC = 0.5 * (ψ[:,:,1] .- ψ[:,:,2])
    ψBT = 0.5 * (ψ[:,:,1] .+ ψ[:,:,2])

    ψBCh = deepcopy(vars.uh[:,:,1])
    ψBTh = deepcopy(vars.uh[:,:,1])

    mul2D!(ψBCh, rfftplan, ψBC)
    mul2D!(ψBTh, rfftplan, ψBT)
    
    U₁, U₂, = view(params.U, :, :, 1), view(params.U, :, :, 2)
    
    S32 = CUDA.@allowscalar f0 * U₁[1,1] / g′

    # streamfunction stuff
    ∂xψBTh = im * grid.kr .* ψBTh
    ∂yψBTh = im * grid.l .* ψBTh
    
    ∂xψBT = deepcopy(vars.u[:,:,1])
    ∂yψBT = deepcopy(vars.u[:,:,1])

    ldiv2D!(∂xψBT, rfftplan, ∂xψBTh)
    ldiv2D!(∂yψBT, rfftplan, ∂yψBTh)

    ∂xψBCh = im * grid.kr .* ψBCh
    ∂yψBCh = im * grid.l .* ψBCh

    ∂xψBC = deepcopy(vars.u[:,:,1])
    ∂yψBC = deepcopy(vars.u[:,:,1])

    ldiv2D!(∂xψBC, rfftplan, ∂xψBCh)
    ldiv2D!(∂yψBC, rfftplan, ∂yψBCh)
    
    # nonlinear terms
    ζ = deepcopy(vars.v)
    ζ[:,:,1] .= irfft(-grid.Krsq .* vars.ψh[:,:,1], grid.ny)
    ζ[:,:,2] .= irfft(-grid.Krsq .* vars.ψh[:,:,2], grid.ny)
    ζ₁, ζ₂ = view(ζ, :, :, 1), view(ζ, :, :, 2)

    ζBCh = - grid.Krsq .* ψBCh
    ζBTh = - grid.Krsq .* ψBTh

    ζBC = deepcopy(vars.u[:,:,1])
    ζBT = deepcopy(vars.u[:,:,1])

    ldiv2D!(ζBC, rfftplan, ζBCh)
    ldiv2D!(ζBT, rfftplan, ζBTh)

    ∂xζBCh = im * grid.kr .* ζBCh
    ∂xζBTh = im * grid.kr .* ζBTh

    ∂yζBCh = im * grid.l .* ζBCh
    ∂yζBTh = im * grid.l .* ζBTh

    ∂xζBC = deepcopy(vars.u[:,:,1])
    ∂xζBT = deepcopy(vars.u[:,:,1])

    ldiv2D!(∂xζBC, rfftplan, ∂xζBCh)
    ldiv2D!(∂xζBT, rfftplan, ∂xζBTh)

    ∂yζBC = deepcopy(vars.u[:,:,1])
    ∂yζBT = deepcopy(vars.u[:,:,1])

    ldiv2D!(∂yζBC, rfftplan, ∂yζBCh)
    ldiv2D!(∂yζBT, rfftplan, ∂yζBTh)

    ##
    ζBT∂xψBTh = deepcopy(vars.uh[:,:,1])
    ζBT∂yψBTh = deepcopy(vars.uh[:,:,1])

    ζBC∂xψBCh = deepcopy(vars.uh[:,:,1])
    ζBC∂yψBCh = deepcopy(vars.uh[:,:,1])

    ζBC∂xψBTh = deepcopy(vars.uh[:,:,1])
    ζBC∂yψBTh = deepcopy(vars.uh[:,:,1])

    ζBT∂xψBCh = deepcopy(vars.uh[:,:,1])
    ζBT∂yψBCh = deepcopy(vars.uh[:,:,1])

    ψBC∂xψBTh = deepcopy(vars.uh[:,:,1])
    ψBC∂yψBTh = deepcopy(vars.uh[:,:,1])

    
    ζBT∂xψBTh = mul2D!(ζBT∂xψBTh, rfftplan, ζBT .* ∂xψBT)
    ζBT∂yψBTh = mul2D!(ζBT∂yψBTh, rfftplan, ζBT .* ∂yψBT)

    ζBC∂xψBCh = mul2D!(ζBC∂xψBCh, rfftplan, ζBC .* ∂xψBC)
    ζBC∂yψBCh = mul2D!(ζBC∂yψBCh, rfftplan, ζBC .* ∂yψBC)

    ζBC∂xψBTh = mul2D!(ζBC∂xψBTh, rfftplan, ζBC .* ∂xψBT)
    ζBC∂yψBTh = mul2D!(ζBC∂yψBTh, rfftplan, ζBC .* ∂yψBT)

    ζBT∂xψBCh = mul2D!(ζBT∂xψBCh, rfftplan, ζBT .* ∂xψBC)
    ζBT∂yψBCh = mul2D!(ζBT∂yψBCh, rfftplan, ζBT .* ∂yψBC)

    ψBC∂xψBTh = mul2D!(ψBC∂xψBTh, rfftplan, ψBC .* ∂xψBT)
    ψBC∂yψBTh = mul2D!(ψBC∂yψBTh, rfftplan, ψBC .* ∂yψBT)

    J_ψBT_ζBT = ∂xψBT .* ∂yζBT .- ∂yψBT .* ∂xζBT

    J_ψBT_ζBTh = deepcopy(vars.uh[:,:,1])
    
    mul2D!(J_ψBT_ζBTh, rfftplan, J_ψBT_ζBT)

    ζ₁h = rfft(ζ₁)
    
    # ζ₂h = rfft(ζ₂)
    
    ∂xζ1h = im * grid.kr .* ζ₁h 
    ∂xζ1 = irfft(∂xζ1h, grid.ny)

    ############################################################################################
    BTEKE = @. 0.5 * (conj(ψBTh) * (grid.kr^2 * ψBTh) + ψBTh * conj(grid.kr^2 * ψBTh))

    BCEKE = @. 0.5 * grid.kr^2 * (conj(ψBCh) * ψBCh + ψBCh * conj(ψBCh))
    BCEAPE = @. 0.5 * 2 * Ld1^-2 * (conj(ψBCh) * ψBCh + ψBCh * conj(ψBCh))

    ############################################################################################
    # Baroclinic conversion term
    CBC = @. S32 * (params.f₀ / params.H[1]) * ( conj(ψBCh) * ∂xψBTh + ψBCh * conj(∂xψBTh)) # these are equivalent!

    # CBC = @. U₁ * Ld1^-2 * (conj(ψBCh) * ∂xψBTh)
    # CBC .+= conj(CBC)
    
    ############################################################################################
    # Linear BC-BT flux (flat-bottom) term
    BC2BT = @. 0.5 * U₁ * (conj(ψBTh) * ∂xζBCh + ψBTh * conj(∂xζBCh))

    ############################################################################################
    # Nonlinear terms in BT budget
    NLBT = zeros(dev, T, (grid.nkr,grid.ny,2)) .+ 0im
    @views NLBT[:,:,1] = @. conj(ψBTh) * J_ψBT_ζBTh + ψBTh * conj(J_ψBT_ζBTh)  # im * (grid.l * ζBT∂xψBTh - grid.kr * ζBT∂yψBTh)
    # NLBT[:,:,1] .+= conj.(NLBT[:,:,1])

    @views NLBT[:,:,2] = @. conj(ψBTh) * im * (grid.l * ζBC∂xψBCh - grid.kr * ζBC∂yψBCh)
    @views NLBT[:,:,2] .+= conj.(NLBT[:,:,2])

    ############################################################################################
    # Nonlinear terms in BC budget
    NLBC = zeros(dev, T, (grid.nkr,grid.ny,3)) .+ 0im

    #NLBC2BT
    @views NLBC[:,:,1] = @. conj(ψBCh) * im * (grid.l * ζBC∂xψBTh - grid.kr * ζBC∂yψBTh)
    @views NLBC[:,:,1] .+= conj.(NLBC[:,:,1])

    #NLBCEKE
    @views NLBC[:,:,2] = @. conj(ψBCh) * im * (grid.l * ζBT∂xψBCh - grid.kr * ζBT∂yψBCh)
    @views NLBC[:,:,2] .+= conj.(NLBC[:,:,2])

    #NLBCEAPE
    @views NLBC[:,:,3] = @. - 2 * Ld1^-2 * conj(ψBCh) * im * (grid.l * ψBC∂xψBTh - grid.kr * ψBC∂yψBTh)
    @views NLBC[:,:,3] .+= conj.(NLBC[:,:,3])

    ############################################################################################
    TopoT = @. -0.5 * (params.f₀ / params.H[2]) * ∂yηb * (conj(ψBTh) * ∂xψBCh + ψBTh * conj(∂xψBCh))

    ############################################################################################
    # BT drag
    DBT = @. - 0.5 * params.μ * (conj(ψBTh) * (ζBCh - ζBTh) + ψBTh * conj(ζBCh - ζBTh))

    ############################################################################################
    # BC drag
    DBC = @. 0.5 * params.μ * (conj(ψBCh) * (ζBCh - ζBTh) + ψBCh * conj(ζBCh - ζBTh) )

    ############################################################################################
    # taking mean at each wavenumber magnitude, i.e., assuming isotropy
    NRGs = hcat(isotropic_mean(BTEKE,grid), isotropic_mean(BCEKE,grid), isotropic_mean(BCEAPE,grid))
    CBCh = isotropic_mean(CBC, grid)
    LF = hcat(isotropic_mean(BC2BT, grid), isotropic_mean(TopoT, grid))
    Drag = hcat(isotropic_mean(DBT, grid), isotropic_mean(DBC, grid))
    NLBTh = hcat(isotropic_mean(NLBT[:,:,1], grid), isotropic_mean(NLBT[:,:,2], grid))  # BC2BT transfer, EKE
    NLBCh = hcat(isotropic_mean(NLBC[:,:,1], grid), isotropic_mean(NLBC[:,:,2], grid), isotropic_mean(NLBC[:,:,3], grid)) # BT2BC, EKE, EAPE

    NRGs = dropdims(NRGs,dims=tuple(findall(size(NRGs).==1)...))
    CBCh = dropdims(CBCh,dims=tuple(findall(size(CBCh).==1)...))
    LF = dropdims(LF,dims=tuple(findall(size(LF).==1)...))
    Drag = dropdims(Drag,dims=tuple(findall(size(Drag).==1)...))
    NLBTh = dropdims(NLBTh,dims=tuple(findall(size(NLBTh).==1)...))
    NLBCh = dropdims(NLBCh,dims=tuple(findall(size(NLBCh).==1)...))

    resid = CBCh .+ sum(Drag, dims=2) .+ sum(NLBTh, dims=2) .+ sum(NLBCh, dims=2)

    ##############################################################################
    ## START: x-space budget
    ##############################################################################
    # Linear modal transfer terms
    LT_x = @. 0.5 * ψBT * U₁ * ∂xζ1   # BC to BT

    # Topo transfer term
    TT_x = - 0.5 * ∂yηb * (params.f₀ / params.H[2]) .* ψBT .* ∂xψBC
    # k-space version: -0.5 * (params.f₀ / params.H[2]) * ∂yηb * (conj(ψBTh) * ∂xψBCh + ψBTh * conj(∂xψBCh))

    # Baroclinic production term
    BC_x = (params.f₀ / params.H[1]) * S32 .* ψBC .* ∂xψBT

    # Nonlinear transfer term
    NL_x = ψBT .* (∂xψBC .* ∂yζBC .- ∂yψBC .* ∂xζBC)  # NL BT <--> BC (EKE)

    # Drag terms
    DBT_x = @. - (μ/2) * ψBT .* (ζBC .- ζBT)     # BT drag term
    DBC_x = @. (μ/2) * ψBC .* (ζBC .- ζBT)      # BC drag term    

    # integrating in y
    LT_x = sum(LT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    TT_x = sum(TT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    BC_x = sum(BC_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    NL_x = sum(NL_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    DBT_x = sum(DBT_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    DBC_x = sum(DBC_x) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    resid_x = BC_x + DBT_x + DBC_x

    # energies
    BTKE_x = 0.5 * sum(∂xψBT.^2 .+ ∂yψBT.^2) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    BCKE_x = 0.5 * sum(∂xψBC.^2 .+ ∂yψBC.^2) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    BCEAPE_x = 2 * Ld1^-2 * sum(ψBC .* ψBC)  * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    ##############################################################################
    ## END: x-space budget
    ##############################################################################

    ## defining KE length scales
    L_BC = sqrt(mean(∂xψBC.^2 .+ ∂yψBC.^2) / mean(ψBC.^2))  # BC
    L_BT = sqrt(mean(∂xψBT.^2 .+ ∂yψBT.^2) / mean(ψBT.^2))  # BT

    GC.gc()
   
    # k-space energies are, BTEKE, BCEKE, EAPE; CBC; Tflat, Ttopo; , DBT, DBC; NLBT2BC, NLBTEKE; NLBC2BT, NLBCEKE, NLBCEAPE resid

  return sqrt(mean(∂xψBT.^2 .+ ∂yψBT.^2)), nrgs_in .+ hcat(NRGs, CBCh, LF, Drag, NLBTh, NLBCh, resid), nrgs_in_x .+ CuArray(vcat(BTKE_x, BCKE_x, BCEAPE_x, LT_x, TT_x, BC_x, NL_x, DBT_x, DBC_x, resid_x)), lengths_in .+ CuArray(vcat(L_BT, L_BC))
end

update_two_layer_kspace_modal_nrgs(prob, ψ, model_params, nrgs_in, nrgs_in_x, lengths_in) = update_two_layer_kspace_modal_nrgs(prob.vars, prob.params, prob.grid, prob.sol, ψ, model_params, nrgs_in, nrgs_in_x, lengths_in)


####################################################################################
## Helpers for k-space budget
####################################################################################

mul2D!(varh, rfftplan, var) = mul!(reshape(varh, (size(varh)...,1)), rfftplan, reshape(var, (size(var)...,1)))[:,:,1]

ldiv2D!(var, rfftplan, varh) = ldiv!(reshape(var, (size(var)...,1)), rfftplan, reshape(varh, (size(varh)...,1)))[:,:,1]


function isotropic_mean(arr_in, grid)
    # arr_in: an nkr X nl array that is output of rfft
    # note that we only want real part of this

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)
    
    dk = 2*pi/grid.Lx; dl = 2*pi/grid.Ly;
    
    dkr = sqrt(dk^2 + dl^2)
    
    wv = @. sqrt(grid.kr^2 + grid.l^2)
    
    iso = zeros(dev, T, (length(grid.kr)))
    
    for i in range(1,length(grid.kr))
        # find 2D index values for a wavenumber magnitude
        if i==length(grid.kr)
            fkr = CUDA.@allowscalar  @. (wv>=grid.kr[i]) & (wv<=grid.kr[i]+dkr)
        else
            fkr = CUDA.@allowscalar  @. (wv>=grid.kr[i]) & (wv<grid.kr[i+1])
        end
        
        if sum(fkr) > 0
            CUDA.@allowscalar iso[i] = mean(real(arr_in[fkr])) # this is average over all combinations of k_x and k_y that are the same, isotropic k
        end
        
    end

    return iso
end

####################################################################################
## x-space layer-wise energies (two-layer for now)
####################################################################################
function update_two_layered_nrg(vars, params, grid, sol, ψ, model_params)

    @unpack_mod_params model_params

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)

    nlayers = 2
    g′ = g * (rho[2] - rho[1]) / rho0[1] # reduced gravity at the interface
    F = (f0^2 / (g′ * sum(params.H)))
    μ = params.μ
    Ld = F^-0.5
    
    ψh = rfft(ψ)

    ∂ψ∂x = irfft(im * grid.kr .* ψh, grid.nx)
    ∂ψ∂y = irfft(im * grid.l  .* ψh, grid.nx)
    
    mod2∇ψ1 = @. ∂ψ∂x[:,:,1]^2 + ∂ψ∂y[:,:,1]^2
    mod2∇ψ2 = @. ∂ψ∂x[:,:,2]^2 + ∂ψ∂y[:,:,2]^2

    # mod2u = @. ∂ψ∂y[:,:,1]^2 + ∂ψ∂y[:,:,2]^2
    # mod2v = @. ∂ψ∂x[:,:,1]^2 + ∂ψ∂x[:,:,2]^2

    APE_d = @. 0.5 * F * (ψ[:,:,1] - ψ[:,:,2])^2
    APE = sum(APE_d) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    
    KE_d = zeros(dev, T, (grid.nx, grid.ny, nlayers))
    @views KE_d[:,:,1] = @. 0.5 * mod2∇ψ1 
    @views KE_d[:,:,2] = @. 0.5 * mod2∇ψ2 
    KE = dropdims(sum(dropdims(sum(KE_d,dims=1),dims=1),dims=1),dims=1) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    # KE_uv = zeros(grid.nx, grid.ny, nlayers)
    # KE_uv[:,:,1] = @. 0.5 * mod2u
    # KE_uv[:,:,2] = @. 0.5 * mod2v 
    # KE_comp = dropdims(sum(dropdims(sum(KE_uv,dims=1),dims=1),dims=1),dims=1) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    
    return CUDA.@allowscalar Array(vcat(KE, APE)) # , nrgs_in .+ vcat(KE, APE)

end

update_two_layered_nrg(prob, ψ, model_params) = update_two_layered_nrg(prob.vars, prob.params, prob.grid, prob.sol, ψ, model_params)

####################################################################################
## calculating the phase shift diagnostics (we use two differnent formulations for now)
####################################################################################

function calc_zonal_phase_shift_12(model_params, psi, grid, vars)
    
    @unpack_mod_params model_params

    k_max_ind = 256

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)

    rfftplan = plan_flows_rfft(A{T, 3}(undef, grid.nx, 1, 1), [1, 2]; flags=FFTW.MEASURE)

    nstep = 20

    nmax = ceil(Int, Ny/nstep)  # 
    
    ph = zeros(nmax,ceil(Int, Nx/2+1))
    ph_wt = zeros(nmax,)

    n=1  # index over all steps

    st = collect(range(1,Ny,step=nstep));

    for k in st   # for each zonal slice where you take an FFT
        # define zonal slice of psi in each layer

        psi1 = psi[:,k,1]
        psi2 = psi[:,k,2]

        psi1h = deepcopy(vars.uh[:,1,1])
        psi2h = deepcopy(vars.uh[:,1,1])
    
        
        mul!(reshape(psi1h, (size(psi1h)...,1,1)), rfftplan, reshape(psi1, (size(psi1)...,1,1)))[:,1,1]
        mul!(reshape(psi2h, (size(psi2h)...,1,1)), rfftplan, reshape(psi2, (size(psi2)...,1,1)))[:,1,1]

        # define weights as a function of wavenumber for this zonal slice
        wts = abs2.(psi1h[1:k_max_ind]) .+ abs2.(psi2h[1:k_max_ind])
        wts = wts ./ sum(wts)
        
        # ph[n,:] = angle.(psi1 ./ psi2)

        ph_wt[n] = sum(wts .* angle.(psi1h[1:k_max_ind] ./ psi2h[1:k_max_ind])) # amplitude-weighted phase at each wavenumber for nth slice

        n+=1
    end

    # ph_out = median(ph[1:n-1,:],dims=1)   # we take median because we want a characteristic phase shift as a function of zonal wavenumber; mean would work too

    ph_wt_out = mean(ph_wt[1:n-1])  # we take mean because this is an array of the summed weighted phase shift where n is the number of zonal slices

    return ph_wt_out
    
end


function calc_iso_phase_shift_12(model_params, psi, grid, vars)

    @unpack_mod_params model_params

    dev = grid.device
    T = eltype(grid)
    A = device_array(dev)

    rfftplan = plan_flows_rfft(A{T, 3}(undef, grid.nx, grid.ny, Nz), [1, 2]; flags=FFTW.MEASURE)

    ψh = deepcopy(vars.uh) # zeros(prob.grid.nkr, prob.grid.ny, Nz) .+ 0im #
    
    mul!(ψh, rfftplan, psi);
    
    wts = isotropic_mean(abs2.(ψh[:,:,1]), grid) .+ isotropic_mean(abs2.(ψh[:,:,2]), grid)
    wts = wts ./ sum(wts)
    
    ph_iso = isotropic_mean(angle.(ψh[:,:,1] ./ ψh[:,:,2]), grid)
    ph_wts = wts .* ph_iso

    return sum(filter(!isnan, ph_wts))
end




####################################################################################
## Diagnosing length scales
####################################################################################







####################################################################################
## 
####################################################################################








####################################################################################
## 
####################################################################################








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

function redef_mu_kappa_topoPV_h0(model_params, mu, kappa, topo_PV, h0_new)
    @unpack_mod_params model_params

    mp_out = mod_params(
    data_dir = data_dir,
    Nz = Nz, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, Ld = Ld,
    H = H,
    rho0 = rho0, rho = rho, strat_str = strat_str,
    shear_str = shear_str, U = U,
    μ = mu , κ = kappa , nν = nν, ν = ν, dyn_nu=dyn_nu,
    eta = eta, topographic_pv_gradient = topo_PV, topo_type = topo_type, h0 = h0_new,
    β = β,
    dt = dt,
    stepper = stepper,
    dev = dev,
    restart_bool = restart_bool,
    restart_yr = restart_yr,
    pre_buoy_restart_file = pre_buoy_restart_file,
    data_dir_pre_buoy = data_dir_pre_buoy,
    ss_yr_max = ss_yr_max,
    yr_increment = yr_increment,
    nsubs = nsubs);

    return mp_out
end

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


function calc_Ld_modes(model_params)

    @unpack_mod_params model_params

    S = calc_stretching_mat(model_params)

    evecs_S = eigvecs(-S)

    # normalizing eigenvectors
    norm_fact = sqrt.(sum(H) ./ (H .* sum( evecs_S .* evecs_S, dims=1)))
    modes = norm_fact .* evecs_S

    return abs.(eigvals(S)).^-0.5, modes
end
