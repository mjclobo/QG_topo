# plotting functions, currently for three layer model
using Measures, PyPlot

function plot_three_layer(tiempo,KE,Cterms,q,grid,kt,h0,plotpath,plotname,ell)

    q1 = transpose(q[:, :, 1])
    q2 = transpose(q[:, :, 2])
    q3 = transpose(q[:, :, 3])
    
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

    pc4=ax4.pcolormesh(grid.x/grid.Lx,grid.y/grid.Ly,q3,cmap=matplotlib.cm.coolwarm,norm=matplotlib.colors.TwoSlopeNorm(0))
    ax4.set_xlabel(L"x/L_x",fontsize=22.)
    ax4.set_ylabel(L"y/L_y",fontsize=22.)
    ax4.set_title(L"q_3",fontsize=26.)
    ax4.tick_params(labelsize=16.)
    cb4 = fig.colorbar(pc4)
    cb4.ax.tick_params(labelsize=16.)

    PyPlot.suptitle(L"k_{topo}= "*string(kt)*L", h_0= "*string(h0),fontsize=30.)

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    PyPlot.savefig(savename)

    PyPlot.close()

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


