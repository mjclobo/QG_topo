# plotting functions, currently for three layer model
using Measures, PyPlot

using PyCall
mpl = pyimport("mpl_toolkits.mplot3d");

function plot_three_layer(tiempo,KE,Cterms,q,v,grid,kt,h0,plotpath,plotname,ell)

    q1 = transpose(q[:, :, 1])
    q2 = transpose(q[:, :, 2])
    q3 = transpose(q[:, :, 3])

    v1 = transpose(v[:, :, 1])
    v2 = transpose(v[:, :, 2])
    v3 = transpose(v[:, :, 3])
    
    fig,ax = PyPlot.subplots(2,2,figsize=(15,10))
    fig.tight_layout(pad=7.0)
    ax1=ax[1]; ax2=ax[2]; ax3=ax[3]; ax4=ax[4];

    ax1.plot(tiempo/3600/24,KE[:,1],linewidth=2.,color="blue",label=L"KE_1 / H_1")
    ax1.plot(tiempo/3600/24,KE[:,2],linewidth=2.,color="orange",label=L"KE_2 / H_2")
    ax1.plot(tiempo/3600/24,KE[:,3],linewidth=2.,color="green",label=L"KE_3 / H_3")
    ax1.tick_params(labelsize=16.)
    ax1.set_xlabel(L"time \quad [days]", fontsize=22.)
    ax1.set_ylabel(L"KE_k / H_k \quad [m s^{-2}]", fontsize=22.)
    # ax1.set_title(L"\psi_1 [\mathrm{norm}]", fontsize = 20.)
    # ax1.set_ylim([t_hovm[1],t_hovm[end]])
    ax1.legend(loc="upper right",fontsize=16.)
  
    p2a = Cterms[1][2:end]
    p2b = Cterms[2][2:end]
    p2c = Cterms[3][2:end] # CL1
    p2T = Cterms[4][2:end] # topography
    p2d = Cterms[5][2:end] # NL1
    p2e = Cterms[6][2:end] # NL2
    p2f = Cterms[7][2:end] # NL3

    ax2.plot(grid.kr[2:end]*grid.Lx/(2*pi),p2a,linewidth=2.,color="blue",label=L"\widehat{C}_{V,3/2}")
    ax2.plot(grid.kr[2:end]*grid.Lx/(2*pi),p2b,linewidth=2.,color="orange",label=L"\widehat{C}_{V,5/2}")
    ax2.plot(grid.kr[2:end]*grid.Lx/(2*pi),p2d,linewidth=2.,color="green",label=L"\widehat{C}_{N,1}")
    ax2.plot([],[],linewidth=2.,color="red",label= L"\widehat{C}_{T}")
    ax2t=ax2.twinx()
    ax2t.plot(grid.kr[2:end]*grid.Lx/(2*pi),p2T,color="red",linewidth=2.,label= L"\widehat{C}_{T} \mathrm{(R)}")
    ax2.set_xlabel(L"k_x L_x /(2 \pi)", fontsize=22.)
    ax2.set_ylabel(L"\mathrm{Spec.} \ \mathrm{energy} \quad [m^3 s^{-3}]",fontsize=22.)
    ax2.legend(loc="upper right",fontsize=16.)
    ax2.tick_params(labelsize=16.)
    ax2t.tick_params(labelsize=16.)

    pc3=ax3.pcolormesh(grid.x/grid.Lx,grid.y/grid.Ly,q1,cmap=matplotlib.cm.coolwarm,norm=matplotlib.colors.TwoSlopeNorm(0))
    ax3.set_xlabel(L"x/L_x",fontsize=22.)
    ax3.set_ylabel(L"y/L_y",fontsize=22.)
    ax3.set_title(L"q_1",fontsize=26.)
    ax3.tick_params(labelsize=16.)
    cb3 = fig.colorbar(pc3)
    cb3.ax.tick_params(labelsize=16.)

    pc4=ax4.pcolormesh(grid.x/grid.Lx,grid.y/grid.Ly,v1,cmap=matplotlib.cm.coolwarm,norm=matplotlib.colors.TwoSlopeNorm(0))
    ax4.set_xlabel(L"x/L_x",fontsize=22.)
    ax4.set_ylabel(L"y/L_y",fontsize=22.)
    ax4.set_title(L"v_1",fontsize=26.)
    ax4.tick_params(labelsize=16.)
    cb4 = fig.colorbar(pc4)
    cb4.ax.tick_params(labelsize=16.)

    PyPlot.suptitle(L"k_{topo}= "*string(kt)*L", h_0= "*string(h0),fontsize=30.)

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    PyPlot.savefig(savename)

    PyPlot.close()

end

function plot_hovm()
    matplotlib[:rcParams]["axes.unicode_minus"]=false

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


end

function plot_growth_rate(k_x,sigma_x,k_emp,sigma_emp,Lx,plotpath)
    mid_int = round(Int,length(k_x)/2+1)
    inv_s_to_day = 3600*24
    p = Plots.plot(k_x[mid_int:end]*Lx/(2*pi),sigma_x[mid_int:end]*inv_s_to_day, linewidth=2.0, xlabel= L"k_x L_x /(2 \pi)", ylabel=L"\sigma \ [ \mathrm{day}^{-1} ]",
                    label="Lin. stab. analysis",xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(750,500),margin=5mm)
    Plots.scatter!([k_emp]*Lx/(2*pi),[sigma_emp]*inv_s_to_day,markersize=6.,color="green",label="Model output")

    local savename = plotpath*"../growth_plot.png"
    Plots.savefig(p,savename)
end

function plot_Qy(H,Qy,plotpath)
    z = [H[1]/2, (H[1]+H[2]/2), (H[1]+H[2]+H[3])/2]
    mq = maximum(abs.(Qy))
    p = Plots.plot(Qy,-z,linewidth=2,ylabel="z [m]",xlabel=L"Q_y",ylimit=[-H[end],0],xlimit=[-mq,mq],
    xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(400,1000),margin=5mm)

    local savename = plotpath*"../Qy.png"
    Plots.savefig(p,savename)
end

function plot_unstable_vert(H,max_eve,psi_emp,plotpath,plotname,ell)
    max_eve_norm = max_eve./maximum(max_eve)
    z = [H[1]/2, (H[1]+H[2]/2), (H[1]+H[2]+H[3])/2]
    p = Plots.plot(max_eve_norm,-z,linewidth=2,ylabel="z [m]",xlabel=L"\psi_\mathrm{norm}",xlimit=[0,1.],ylimit=[-H[end],0],
    # yaxis=:log,
    label=L"\psi_\mathrm{max} \ (L.S.A.)",xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(400,1000),margin=5mm)
    Plots.plot!(psi_emp,-z,linewidth=2,label=L"\psi_\mathrm{max} \ (model)")

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    Plots.savefig(p,savename)

end

function plot_layerwise_spectra(k,specs,plotpath,plotname,ell)

    p = Plots.plot(k[2:end], specs[2:end,:],label=[L"\widehat{\psi}_{1}" L"\widehat{\psi}_{2}" L"\widehat{\psi}_{3}"],
    xlabel= L"k_x L_x /(2 \pi)" ,ylabel=L"\psi \ \mathrm{ PSD} \ [m^4 s^{-2} m]",linewidth=2., linecolor=[:green :red :blue], 
    xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(750,500),margin=5mm)

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    Plots.savefig(p,savename)

end

# calculate phase speeds
function calc_phase_speeds(psi_ot,t_hovm,qy1,U,Lx,Nx)
    psi1_ot = psi_ot[1]; psi2_ot = psi_ot[2]; psi3_ot = psi_ot[3]; 
    x = collect(range(0.0,Lx,Nx))
    dpsi1_dt = diff(psi1_ot,dims=2) ./ diff(t_hovm)'
    dpsi1_dx = diff(psi1_ot,dims=1)' ./ diff(x)'
    dpsi2_dt = diff(psi2_ot,dims=2) ./ diff(t_hovm)'
    dpsi2_dx = diff(psi2_ot,dims=1)' ./ diff(x)'
    dpsi3_dt = diff(psi3_ot,dims=2) ./ diff(t_hovm)'
    dpsi3_dx = diff(psi3_ot,dims=1)' ./ diff(x)'

    cr1_dopp = median(abs.(dpsi1_dt[1:end-1,:]./dpsi1_dx[1:end-1,:]'))
    cr2_dopp = median(abs.(dpsi2_dt[1:end-1,:]./dpsi2_dx[1:end-1,:]'))
    cr3_dopp = median(abs.(dpsi3_dt[1:end-1,:]./dpsi3_dx[1:end-1,:]'))

    # positive PV gradient means westward phase speed
    # do we need to divide by 2pi? I don't think so
    # cr1 = qy1[1] > 0.0 ? cr1_dopp - U[1] : U[1] + cr1_dopp
    # cr2 = qy1[2] > 0.0 ? cr2_dopp - U[2] : U[2] + cr2_dopp
    # cr3 = qy1[3] > 0.0 ? cr3_dopp - U[3] : U[3] + cr3_dopp
    cr1 = cr1_dopp - U[1]
    cr2 = cr2_dopp - U[2]
    cr3 = cr3_dopp - U[3]

    # cr1_dopp = cr1; cr2_dopp = cr2; cr3_dopp = cr3

    return [cr1, cr2, cr3], [cr1_dopp, cr2_dopp, cr3_dopp]
end

# calculate CSP criterion, outputting terms at various steps
function calc_csp_crit(cr,qy,U,psi,Nz)
    # assumes constancy in meridional direction
    # phase speeds are Doppler-shifted
    terms = zeros(Nz)
    for i=range(1,Nz)
        terms[i] = H[i] * qy[i] * psi[i]^2 / ((U[i]-cr[i])^2)
    end
    csp_crit = sum(terms)
    return csp_crit, terms
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function plot_box(psi1_full,psi2_full,psi3_full,Lx,Nx,h,plotpath,plotname,ell)

    psi_west = [psi1_full[1,:]'; psi1_full[1,:]'; psi2_full[1,:]'; psi2_full[1,:]'; psi3_full[1,:]'; psi3_full[1,:]']'

    psi_south = [psi1_full[:,1]'; psi1_full[:,1]'; psi2_full[:,1]'; psi2_full[:,1]'; psi3_full[:,1]'; psi3_full[:,1]']'

    psi_all = [psi1_full; psi2_full; psi3_full];

    cmap = PyPlot.cm.bwr
    norm_t = matplotlib.colors.TwoSlopeNorm(0,vmin=minimum(psi_all),vmax=maximum(psi_all))
    norm_w = matplotlib.colors.TwoSlopeNorm(0,vmin=minimum(psi_all),vmax=maximum(psi_all))
    norm_s = matplotlib.colors.TwoSlopeNorm(0,vmin=minimum(psi_all),vmax=maximum(psi_all))


    colors_psi_t = cmap(norm_t(psi1_full));
    colors_psi_w = cmap(norm_w(psi_west));
    colors_psi_s = cmap(norm_s(psi_south));

    x = y = collect(range(-Lx/2,Lx/2,Nx))

    eps = 0.01
    z = [0.0,-500. + eps,-500. - eps, -1500. + eps, -1500. - eps, -4000.] .* 10

    ax1 = PyPlot.figure(figsize=(10,4)).add_subplot(projection="3d")

    # west face
    yw,zw = meshgrid(y,z)

    ax1.plot_surface(-Lx/2*ones(size(yw)),yw,zw,facecolors=colors_psi_w,shade=false)

    # top
    xt,yt = meshgrid(x,y)
    ax1.plot_surface(xt,yt,10 .+ zeros(size(xt)),facecolors=colors_psi_t,shade=false)

    # south face
    xs,zs = meshgrid(x,z)

    ax1.plot_surface(xs,-Lx/2*ones(size(xs)),zs,facecolors=colors_psi_s,shade=false)

    # sloping bottom
    codes_sq = [matplotlib.path.Path.MOVETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.CLOSEPOLY];
    codes = [matplotlib.path.Path.MOVETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.LINETO, matplotlib.path.Path.CLOSEPOLY];
    if h>0
        verts = [[-Lx/2,-40000],[Lx/2,-40000],[Lx/2,-40000+10*h*Lx],[-Lx/2,-40000.]]
        pp = matplotlib.patches.PathPatch(matplotlib.path.Path(verts,codes),ec="none",color="gray")
        ax1.add_patch(pp)
        mpl.art3d.pathpatch_2d_to_3d(pp, z=-Lx/2, zdir="x")
    elseif h<0
        verts = [[Lx/2,-40000],[-Lx/2,-40000-10*h*Lx],[-Lx/2,-40000],[Lx/2,-40000]]
        pp = matplotlib.patches.PathPatch(matplotlib.path.Path(verts,codes),ec="none",color="gray")
        ax1.add_patch(pp)
        mpl.art3d.pathpatch_2d_to_3d(pp, z=-Lx/2, zdir="x")

        verts_sq = [[-Lx/2,-40000-10*h*Lx],[-Lx/2,-40000],[Lx/2,-40000.],[Lx/2,-40000-10*h*Lx],[-Lx/2,-40000-10*h*Lx]]
        pp = matplotlib.patches.PathPatch(matplotlib.path.Path(verts_sq,codes_sq),ec="none",color="gray")
        ax1.add_patch(pp)
        mpl.art3d.pathpatch_2d_to_3d(pp, z=-Lx/2, zdir="y")

    else
        a=1.
    end

    # verts = [[0,0],[1,0],[1,1],[0,0]]
    # p2 = matplotlib.patches.Wedge((0,0), 1,0,20, alpha=1.0)

    ax1.view_init(elev=20., azim=250);

    PyPlot.axis("off");

    ax1.set_ylim(-Lx/2,Lx/2);
    
    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    PyPlot.savefig(savename,bbox_inches="tight")

    PyPlot.close()

end

## function to calculate normal modes of system
function vert_modes(rho,h0,g,f0,Lx,Nx)
    Nz = length(rho)

    Ny = Nx
    Ly = Lx

    # define wavenumbers
    k_x = reshape(fftfreq(Nx, 2π/Lx*Nx),(1,Nx))
    k_y = reshape(fftfreq(Ny, 2π/Ly*Ny),(1,Ny))
    k2  = k_x.^2 + k_y.^2 

    N_x = Nx
    N_xr = round(Int,Nx/2)
    k_xr = k_x[1:N_xr]
    k_xr2 = k_xr.^2

    # define stretching matrix
    rigid_lid = true
    S = LinStab.calc_stretching_mat(Nz,rho,f0,H,rho[1],g,rigid_lid,h0)
    
    # def rad
    evals_S = eigvals(-S); evecs_S = eigvecs(-S)
    sort_ind = sortperm(abs.(evals_S))

    # T2 = eltype(eigvals)
    r_d = Complex.(zeros(Nz,1))
    r_d[1] = sqrt(g*sum(H))/f0
    r_d[2:end] = @. sqrt(Complex(evals_S[sort_ind[2:end]]))^-1

    # normalizing eigenvectors
    norm_fact = sqrt.(sum(H) ./ (H .* sum( evecs_S .* evecs_S, dims=1)))
    evecs_S = norm_fact .* evecs_S
    
    return evecs_S, r_d, k_xr2, k_xr, S
end
    
## function to project slice of field onto vertical modes
function project_onto_modes(modes,field)
    # slice of field is Nx x Nz center slice in domain
    # modes is Nz x Nz
    fieldh = rfft(field,2)
    proj = modes \ fieldh
    return proj
end

## function to calculate zonal modal KE spec
function modal_nrg_spec(kx2,S,A)
    # k2 is Nx x 1
    # S is Nz x Nz
    # A is Nx x Nz
    Nx = length(kx2)
    Nz = size(S)[1]
    
    Id = Matrix{Float64}(I, Nz,Nz) # identity matrix
    
    # inversion matrix is Nz x Nz x Nx x Ny
    spec = zeros(Nx,Nz)
    for i=range(1,Nx)
        for j=range(1,Nz)
            
#             M = S .- kx2[i] * Id
            M = Nx^2
    
            spec[i,j] = kx2[i] * ((abs(A[j,i])^2) / (M^2))
            
        end
    end
    
    return spec
end

## integrate modal spectrum to find total energy in each mode
function spec_integration(k_xr,spec)
    # spec is Nx x Nz
    # k_xr is Nx/2
    
    dk = k_xr[2] - k_xr[1]
    
    modal_amp = zeros(Nz,)
    for i=range(1,Nz)
        modal_amp[i] = dk*sum(spec[:,i])
    end
    
    return modal_amp
end


