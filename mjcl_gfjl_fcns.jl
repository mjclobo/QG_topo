# This is a module of functions used for running QG models using GeophysicalFlows.jl.
# 

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
    restart_yr::Int = 0
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

    @unpack_diag_bools diags

    @unpack_mod_params model_params

    sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid

    startwalltime = time()

    global j = 0
    global t_yrly = nothing
    global yr_cnt = restart_yr
    global budget_counter = 0

    nterms_two_layer_modal_kspace = 14  # including residual
    global two_layer_kspace_modal_nrgs = zeros(grid.nkr, nterms_two_layer_modal_kspace)
    global two_layer_xspace_layer_nrgs = zeros(3)
    global two_layer_vBT_scale = 0.

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

            if psi_out_bool==true
                if isnothing(t_yrly)
                    global psi_ot = deepcopy(vars.ψ);
                else
                    global psi_ot = cat(psi_ot, vars.ψ, dims=4)
                end
            end

            if isnothing(t_yrly)
                global t_yrly = Array([clock.t])
            else
                push!(t_yrly,clock.t)
            end

            if two_layer_kspace_modal_nrg_budget_bool==true
                
                global two_layer_kspace_modal_nrgs = update_two_layer_kspace_modal_nrgs(prob, vars.ψ, model_params, two_layer_kspace_modal_nrgs)
                global two_layer_xspace_layer_nrgs = update_two_layered_nrg(prob, vars.ψ, model_params, two_layer_xspace_layer_nrgs) 
                global two_layer_vBT_scale += sqrt(mean((sum(vars.u, dims=2) ./ 2).^2))
                global budget_counter +=1

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
                            "psi_ot" => Array(psi_ot), "two_layer_kspace_modal_nrg_budget" => two_layer_kspace_modal_nrgs)
                    elseif psi_out_bool==false && two_layer_kspace_modal_nrg_budget_bool==true
                        jld_data = Dict("two_layer_kspace_modal_nrg_budget" => two_layer_kspace_modal_nrgs ./ budget_counter,
                            "two_layer_xspace_layer_nrgs" => two_layer_xspace_layer_nrgs ./ budget_counter,
                            "two_layer_vBT_scale" => two_layer_vBT_scale ./ budget_counter)
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


# function set_output_dict(diags)

# end

####################################################################################
## DIAGS: x-space modal budget
####################################################################################






####################################################################################
## DIAGS: k-space budget
####################################################################################

function update_two_layer_kspace_modal_nrgs(vars, params, grid_jl, sol, ψ, model_params, nrgs_in)
    # energies are: BTEKE, BCEKE, EAPE; CBC, DBC, DBT; Tflat, Ttopo; NLBCEAPE, NLBCEKE, NLBC2BT; NLBTEKE, NLBT2BC; resid
    # here we do not define average, just add up the budget...averaging comes later

    @unpack_mod_params model_params

    # parameters
    nlayers = 2
    δ = [params.H[1]/sum(params.H), params.H[2]/sum(params.H)]
    g′ = g * (rho[2] -   rho[1]) / rho0 # reduced gravity at each interface
    Ld1 = (f0^2 / (g′ * params.H[1]))^-0.5
    ∂yηb = params.topographic_pv_gradient[2] / (f0 / params.H[end])
    
    # assigning basic variables
    ψBC = 0.5 * (ψ[:,:,1] .- ψ[:,:,2])
    ψBT = 0.5 * (ψ[:,:,1] .+ ψ[:,:,2])
    
    ψBCh = rfft(ψBC)
    ψBTh = rfft(ψBT)
    
    U₁, U₂, = view(params.U, :, :, 1), view(params.U, :, :, 2)
    
    S32 = CUDA.@allowscalar f0 * U₁[1,1] / g′

    # streamfunction stuff
    ∂xψBTh = im * grid_jl.kr .* ψBTh
    ∂yψBTh = im * grid_jl.l .* ψBTh
    
    ∂xψBT = irfft(∂xψBTh, grid_jl.ny)
    ∂yψBT = irfft(∂yψBTh, grid_jl.ny)

    ∂xψBCh = im * grid_jl.kr .* ψBCh
    ∂yψBCh = im * grid_jl.l .* ψBCh
    
    ∂xψBC = irfft(∂xψBCh, grid_jl.ny)
    ∂yψBC = irfft(∂yψBCh, grid_jl.ny)
    
    # nonlinear terms
    # ζ = deepcopy(vars.v)
    # ζ[:,:,1] .= irfft(-grid_jl.Krsq .* vars.ψh[:,:,1], grid_jl.ny)
    # ζ[:,:,2] .= irfft(-grid_jl.Krsq .* vars.ψh[:,:,2], grid_jl.ny)
    # ζ₁, ζ₂ = view(ζ, :, :, 1), view(ζ, :, :, 2)

    ζBCh = - grid_jl.Krsq .* ψBCh
    ζBTh = - grid_jl.Krsq .* ψBTh

    ζBC = irfft(ζBCh,grid_jl.ny)
    ζBT = irfft(ζBTh,grid_jl.ny)

    ∂xζBCh = im * grid_jl.kr .* ζBCh
    ∂xζBTh = im * grid_jl.kr .* ζBTh

    ∂yζBCh = im * grid_jl.l .* ζBCh
    ∂yζBTh = im * grid_jl.l .* ζBTh

    ∂xζBC = irfft(∂xζBCh, grid_jl.ny)
    ∂xζBT = irfft(∂xζBTh, grid_jl.ny)

    ∂yζBC = irfft(∂yζBCh, grid_jl.ny)
    ∂yζBT = irfft(∂yζBTh, grid_jl.ny)
    
    ζBT∂xψBTh = rfft(ζBT .* ∂xψBT)
    ζBT∂yψBTh = rfft(ζBT .* ∂yψBT)

    ζBC∂xψBCh = rfft(ζBC .* ∂xψBC)
    ζBC∂yψBCh = rfft(ζBC .* ∂yψBC)

    ζBC∂xψBTh = rfft(ζBC .* ∂xψBT)
    ζBC∂yψBTh = rfft(ζBC .* ∂yψBT)

    ζBT∂xψBCh = rfft(ζBT .* ∂xψBC)
    ζBT∂yψBCh = rfft(ζBT .* ∂yψBC)

    ψBC∂xψBT = rfft(ψBC .* ∂xψBT)
    ψBC∂yψBT = rfft(ψBC .* ∂yψBT)

    J_ψBT_ζBT = ∂xψBT .* ∂yζBT .- ∂yψBT .* ∂xζBT

    J_ψBT_ζBTh = rfft(J_ψBT_ζBT)

    # ζ₁h = rfft(ζ₁)
    
    # ζ₂h = rfft(ζ₂)
    
    # ∂xζ1h = im * grid_jl.kr .* ζ₁h 
    # ∂xζ1 = irfft(∂xζ1h, grid_jl.ny)

    ############################################################################################
    BTEKE = @. 0.5 * (conj(ψBTh) * (grid_jl.kr^2 * ψBTh) + ψBTh * conj(grid_jl.kr^2 * ψBTh))

    BCEKE = @. 0.5 * grid_jl.kr^2 * (conj(ψBCh) * ψBCh + ψBCh * conj(ψBCh))
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
    NLBT = zeros(grid_jl.nkr,grid_jl.ny,2) .+ 0im
    @views NLBT[:,:,1] = @. conj(ψBTh) * J_ψBT_ζBTh + ψBTh * conj(J_ψBT_ζBTh)  # im * (grid_jl.l * ζBT∂xψBTh - grid_jl.kr * ζBT∂yψBTh)
    # NLBT[:,:,1] .+= conj.(NLBT[:,:,1])

    @views NLBT[:,:,2] = @. conj(ψBTh) * im * (grid_jl.l * ζBC∂xψBCh - grid_jl.kr * ζBC∂yψBCh)
    @views NLBT[:,:,2] .+= conj.(NLBT[:,:,2])

    ############################################################################################
    # Nonlinear terms in BC budget
    NLBC = zeros(grid_jl.nkr,grid_jl.ny,3) .+ 0im
    @views NLBC[:,:,1] = @. conj(ψBCh) * im * (grid_jl.l * ζBC∂xψBTh - grid_jl.kr * ζBC∂yψBTh)
    @views NLBC[:,:,1] .+= conj.(NLBC[:,:,1])

    @views NLBC[:,:,2] = @. conj(ψBCh) * im * (grid_jl.l * ζBT∂xψBCh - grid_jl.kr * ζBT∂yψBCh)
    @views NLBC[:,:,2] .+= conj.(NLBC[:,:,2])

    @views NLBC[:,:,3] = @. - 2 * Ld1^-2 * conj(ψBCh) * im * (grid_jl.l * ψBC∂xψBT - grid_jl.kr * ψBC∂yψBT)
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
    NRGs = hcat(isotropic_mean(BTEKE,grid_jl), isotropic_mean(BCEKE,grid_jl), isotropic_mean(BCEAPE,grid_jl))
    CBCh = isotropic_mean(CBC, grid_jl)
    LF = hcat(isotropic_mean(BC2BT, grid_jl), isotropic_mean(TopoT, grid_jl))
    Drag = hcat(isotropic_mean(DBT, grid_jl), isotropic_mean(DBC, grid_jl))
    NLBTh = hcat(isotropic_mean(NLBT[:,:,1], grid_jl), isotropic_mean(NLBT[:,:,2], grid_jl))  # BC2BT transfer, EKE
    NLBCh = hcat(isotropic_mean(NLBC[:,:,1], grid_jl), isotropic_mean(NLBC[:,:,2], grid_jl), isotropic_mean(NLBC[:,:,3], grid_jl)) # BT2BC, EKE, EAPE

    NRGs = dropdims(NRGs,dims=tuple(findall(size(NRGs).==1)...))
    CBCh = dropdims(CBCh,dims=tuple(findall(size(CBCh).==1)...))
    LF = dropdims(LF,dims=tuple(findall(size(LF).==1)...))
    Drag = dropdims(Drag,dims=tuple(findall(size(Drag).==1)...))
    NLBTh = dropdims(NLBTh,dims=tuple(findall(size(NLBTh).==1)...))
    NLBCh = dropdims(NLBCh,dims=tuple(findall(size(NLBCh).==1)...))

    resid = CBCh .+ sum(Drag, dims=2) .+ sum(NLBTh, dims=2) .+ sum(NLBCh, dims=2)
    
    GC.gc()
   
    # energies are, BTEKE, BCEKE< EAPE; CBC, DBC, DBT; Tflat, Ttopo; NLBCEAPE, NLBCEKE, NLBC2BT; NLBTEKE, NLBT2BC; resid

  return nrgs_in .+ hcat(NRGs, CBCh, LF, Drag, NLBTh, NLBCh, resid)
end

update_two_layer_kspace_modal_nrgs(prob, ψ, model_params, nrgs_in) = update_two_layer_kspace_modal_nrgs(prob.vars, prob.params, prob.grid, prob.sol, ψ, model_params, nrgs_in)


####################################################################################
## Helper for k-space budget
####################################################################################

function isotropic_mean(arr_in, grid)
    # arr_in: an nkr X nl array that is output of rfft
    # note that we only want real part of this
    
    dk = 2*pi/grid.Lx; dl = 2*pi/grid.Ly;
    
    dkr = sqrt(dk^2 + dl^2)
    
    wv = @. sqrt(grid.kr^2 + grid.l^2)
    
    iso = zeros(length(grid.kr))
    
    for i in range(1,length(grid.kr))
        # find 2D index values for a wavenumber magnitude
        if i==length(grid.kr)
            fkr = @. (wv>=grid.kr[i]) & (wv<=grid.kr[i]+dkr)
        else
            fkr = @. (wv>=grid.kr[i]) & (wv<grid.kr[i+1])
        end
        
        if sum(fkr) > 0
            iso[i] = mean(real(arr_in[fkr])) # this is average over all combinations of k_x and k_y that are the same, isotropic k
        end
        
    end

    return iso
end

####################################################################################
## x-space layer-wise energies (two-layer for now)
####################################################################################
function update_two_layered_nrg(vars, params, grid, sol, ψ, model_params, nrgs_in)

    @unpack_mod_params model_params

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
    
    KE_d = zeros(grid.nx, grid.ny, nlayers)
    KE_d[:,:,1] = @. 0.5 * mod2∇ψ1 
    KE_d[:,:,2] = @. 0.5 * mod2∇ψ2 
    KE = dropdims(sum(dropdims(sum(KE_d,dims=1),dims=1),dims=1),dims=1) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1

    # KE_uv = zeros(grid.nx, grid.ny, nlayers)
    # KE_uv[:,:,1] = @. 0.5 * mod2u
    # KE_uv[:,:,2] = @. 0.5 * mod2v 
    # KE_comp = dropdims(sum(dropdims(sum(KE_uv,dims=1),dims=1),dims=1),dims=1) * grid.dx * grid.dy * grid.Lx^-1 * grid.Ly^-1
    
    return nrgs_in .+ hcat(KE[1], KE[2], APE)

end

update_two_layered_nrg(prob, ψ, model_params, nrgs_in) = update_two_layered_nrg(prob.vars, prob.params, prob.grid, prob.sol, ψ, model_params, nrgs_in)

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
    μ = mu , κ = kappa , nν = nν, ν = ν,
    eta = eta, topographic_pv_gradient = topo_PV, topo_type = topo_type, h0 = h0_new,
    β = β,
    dt = dt,
    stepper = stepper,
    dev = dev,
    restart_bool = restart_bool,
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
