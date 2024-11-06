# This is the kernel for running QG models using GeophysicalFlows.jl.

global H_bound = (H[1:end-1] .+ H[2:end]) ./ 2

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

topographic_pv_gradient = (0., 0.)
aliased_fraction=1/3; T=Float64;
grid_topo = TwoDGrid(dev; nx=Nx, Lx, ny=Ny, Ly, aliased_fraction, T)
if topo_type=="eggshell"
    eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"eggshell")
elseif topo_type=="sinusoid"
    eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sinusoid")
elseif topo_type=="sin_sin"
    eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"sin_sin")
elseif topo_type=="y_slope"
    # eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"y_slope")
    eta = nothing
    topographic_pv_gradient = (0., h0*(f0/H[end]))
elseif topo_type=="rand_slope"
    eta = topographicPV(grid_topo,h0[2],kt,Lx,Ly,f0,H,"rand")
    topographic_pv_gradient = (0., h0[1]*(f0/H[end]))
elseif topo_type=="rand_flat"
    eta = topographicPV(grid_topo,h0,kt,Lx,Ly,f0,H,"rand")
    topographic_pv_gradient = (0., 0.)
else
    eta = nothing
end

# define the model problem
# prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, g, H, ρ, U, eta=eta, nν, ν,
# μ, β, dt, stepper, linear, aliased_fraction=1/3)

b = (g/rho0)*(rho0 .- rho)

# prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, H, b, U, nν, ν, eta, topographic_pv_gradient,
# μ, β, dt, stepper, linear, aliased_fraction=1/3)

if ν==0.
    af = 0
else
    af=1/3
end

prob = MultiLayerQG.Problem(nlayers, dev; nx=n, Lx=L, f₀, H, g, ρ, U, nν, ν, eta, topographic_pv_gradient,
μ, β, dt, stepper, linear, aliased_fraction=af, drag_bool)

sol, clock, params, vars, grid = prob.sol, prob.clock, prob.params, prob.vars, prob.grid
x, y = grid.x, grid.y

# setting initial conditions
seed!(1234) # reset of the random number generator for reproducibility
q₀  = q0_mag * device_array(dev)(randn((grid.nx, grid.ny, nlayers)))
q₀h = prob.timestepper.filter .* rfft(q₀, (1, 2)) # apply rfft  only in dims=1, 2
q₀  = irfft(q₀h, grid.nx, (1, 2))                 # apply irfft only in dims=1, 2

MultiLayerQG.set_q!(prob, q₀)

# output dirs
filepath = "."

if topo_type=="y_slope"
    plotpath_main = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/main/"
    plotpath_psi  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi/"
    plotpath_psi_vert  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0*Lx,digits=9))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
elseif topo_type=="sin_sin"
    plotpath_main = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/main/"
    plotpath_psi  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/psi/"
    plotpath_psi_vert  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
elseif topo_type=="rand_slope"
    plotpath_main = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/main/"
    plotpath_psi  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/psi/"
    plotpath_psi_vert  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(round(h0[1]*Lx,digits=9))*"_kt"* string(Int(kt[1])) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
elseif linear
    plotpath_main = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/main/"
    plotpath_psi  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi/"
    plotpath_psi_vert  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_linear_res" * string(Int(Nx)) *"/psi_vert/"
else
    plotpath_main = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/main/"
    plotpath_psi  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi/"
    plotpath_psi_vert  = "./figs/plots_"*string(nlayers)*"layer_h0"* string(Int(h0))*"_kt"* string(Int(kt)) *"_res" * string(Int(Nx)) *"/psi_vert/"
end
plotname = "snapshots"
include("./plotting_functions.jl")

# file management
if !isdir(plotpath_main); mkpath(plotpath_main); end
if !isdir(plotpath_psi); mkpath(plotpath_psi); end
if !isdir(plotpath_psi_vert); mkpath(plotpath_psi_vert); end

# trying to run the model now
startwalltime = time()

global ell, j = 1, 0
global t_yrly = nothing
global yr_cnt = 0
global ss_yr = false
global ss_yr_cnt = 0

while ss_yr_cnt < ss_yr_max
    global ell, j 

    ##########################
    stepforward!(prob)
    MultiLayerQG.updatevars!(prob)

    if j % nsubs == 0
        local E = MultiLayerQG.energies(prob)

        if isnothing(t_yrly)
            global psi_ot = vars.ψ;

            global t_yrly = Array([clock.t])

            # defining initial diagnostics
            E = MultiLayerQG.energies(prob)

            # variables for plotting, that will be pushed to
            global KE = E[1] ./ H

            global PE = E[2] ./ H_bound

        else

            global psi_ot = cat(psi_ot, vars.ψ, dims=4)
            
            # defining initial diagnostics
            E = MultiLayerQG.energies(prob)

            # push to time
            push!(t_yrly,clock.t)

            # push to layer-wise KE
            global KE = cat(KE, E[1] ./ H, dims=2)

            # push to interface potential energies
            global PE = cat(PE, E[2] ./ H_bound, dims=2)

            # plotting stuff
            global plot_model
            if plot_model==true

            end

            # increase counter
            global ell+=1

        end
    
        # reading out stats
        cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
    
        log = @sprintf("step: %04d, t: %.1f, cfl: %.2f, KE_1: %.3e, KE_N: %.3e, PE: %.3e, walltime: %.2f min",
                        clock.step, clock.t, cfl, E[1][1], E[1][end], E[2][1], (time()-startwalltime)/60)
        
        println(log)

        if E[1][1]==NaN
            global ss_yr_cnt = ss_yr_max
        end

        # save output and reset params every year
        if ((t_yrly[end] - yr_cnt*365*86400) > 0)
            global yr_cnt += 1
            global ell = 0
            # check to see if model has reached s.s.
            if ss_yr==true
                global ss_yr_cnt += 1
            end

            if KE[1,end]/KE[1,1] < KE_thresh && ss_yr==false && yr_cnt > 10
                global ss_yr = true
                global ss_yr_cnt = 1
            end

            if isnan(KE[1,end])
                global ss_yr_cnt = ss_yr_max
            end

            # saving yearly output
            println("Saving annual data for year: "*string(yr_cnt))

            if drag_bool==true
                drag_str="_quad_drag_"
                mu_str = @sprintf "%.2E" μ * Ld
            else
                drag_str="_lin_drag_"
                mu_str = @sprintf "%.2E" μ * 2 * Ld / U[1]
            end

            nu_str = @sprintf "%.2E" ν / ((Us[1]/2) * (Lx/2/pi)^7)

            if topo_type=="y_slope"
                jld_name = data_dir*"/threelayer_h0"* string(round(h0/S32,digits=3))* "_beta" * string(round(β * 2 * Ld^2 / U[1],digits=3)) * "_U" * string(round(U[1],digits=4)) * "_rho"* string(round(ρ[1],digits=6)) * drag_str * "mu" * mu_str * "_nu" * nu_str * "_Hr" * string(round(H[1]/H[2],sigdigits=1)) * "_res" * string(Int(Nx)) * "_yr" * string(yr_cnt) *  ".jld"
            elseif topo_type=="rand_slope"
                jld_name = data_dir*"/threelayer_h0"* string(round(h0[1]*Lx,digits=9))* "_hrms"* string(round(h0[2]))*"_kt" * string(round(kt)) * "_U" * string(round(U[1],digits=6)) * "_rho"* string(round(ρ[1],digits=6)) * lin_str * "mu" * string(round((μ^-1)/86400)) * "Hr" * string(round(H[1]/H[2],sigdigits=1))  * "res" * string(Int(Nx)) * "_yr"*string(yr_cnt)*  ".jld"
            elseif topo_type=="sin_sin"
                jld_name = data_dir*"/threelayer_h0"* string(Int(h0[1]))* "_U" * string(round(U[1],digits=6)) * "_rho"* string(round(ρ[1],digits=6)) *lin_str * "mu" * string(round((μ^-1)/86400)) * "Hr" * string(round(H[1]/H[2],sigdigits=1))  * "_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld"  
            else
                jld_name = data_dir*"/threelayer_h0"* string(Int(h0))* "_U" * string(round(U[1],digits=6)) * "_rho"* string(round(ρ[1],digits=6)) * lin_str * "mu" * string(round((μ^-1)/86400)) * "Hr" * string(round(H[1]/H[2],sigdigits=1))  * "_res" * string(Int(Nx)) *"_yr"*string(yr_cnt)*  ".jld"
            end

            println("Saving output data to JLD to: "*jld_name)

            jld_data = Dict("t" => t_yrly, "KE" => KE, "Nz" => nlayers,
                            "L" => L, "H" => H, "rho" => rho, "U" => U,
                            "dt" => dt, "beta" => β, "mu" => μ,
                            "psi_ot" => Array(psi_ot),
                            "nu" => ν, "n_nu" => nν, "Ld" => Ld,
                            "Qy" => Array(params.Qy), "eta" => eta,
                            "cfl_set" => cfl_glob, "PE" => PE)
        
            jldsave(jld_name; jld_data)

            global t_yrly = nothing
        end

        GC.gc()

    end

    global j+=1
    ###############################

end
