# plotting functions, currently for three layer model
using Measures

function plot_three_layer(tiempo,KE,Cterms,q,grid,kt,h0,plotpath,plotname,ell)

    q1 = transpose(q[:, :, 1])
    q2 = transpose(q[:, :, 2])
    q3 = transpose(q[:, :, 3])
    
    p1 = plot(tiempo, KE,label=["layer 1" "layer 2" "layer 3"],xlabel="time",ylabel=L"KE [m s^{-2}]")
  
    p2a = Cterms[1][2:end]

    p2b = Cterms[2][2:end]

    p2c = Cterms[3][2:end]

    p2 = plot(grid.kr[2:end]*grid.Lx/(2*pi), [p2a p2b p2c],label=[L"\widehat{C}_{V,3/2}" L"\widehat{C}_{V,5/2}" L"\widehat{C}_{L,1}"],
          xaxis=:log, xlabel= L"k_x L_x /(2 \pi)" ,ylabel=L"\mathrm{Spectral thickness flux }[m^3 s^{-3}]")

    p3 = heatmap(grid.x/grid.Lx,grid.y/grid.Ly,q1,xlabel=L"x/L_x ",ylabel=L"y/L_y",xflip=true,title=L"q_1",c=:balance,colorbar_exponentformat="power")
    p4 = heatmap(grid.x/grid.Lx,grid.y/grid.Ly,q3,xlabel=L"x/L_y",ylabel=L"y/L_y",xflip=true,title=L"q_3",c=:balance,colorbar_exponentformat="power")

    p = plot(p1,p2,p3,p4,layout=(2,2),plot_title=L"k_{topo}= "*string(kt)*L", h_0= "*string(h0)*L", Bu_{1,2} = "*type,size=(1200,750))

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    savefig(p,savename)

end


function plot_growth_rate(k_x,sigma_x,k_emp,sigma_emp,Lx,plotpath)
    mid_int = round(Int,length(k_x)/2+1)
    inv_s_to_day = 3600*24
    p = Plots.plot(k_x[mid_int:end]*Lx/(2*pi),sigma_x[mid_int:end]*inv_s_to_day, linewidth=2.0, xlabel= L"k_x L_x /(2 \pi)", ylabel=L"\sigma \ [ \mathrm{day}^{-1} ]",
                    label="Lin. stab. analysis",xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(750,500),margin=5mm)
    Plots.scatter!([k_emp]*Lx/(2*pi),[sigma_emp]*inv_s_to_day,markersize=6.,color="green",label="Model output")

    local savename = plotpath*"../growth_plot.png"
    savefig(p,savename)
end

function plot_Qy(H,Qy,plotpath)
    z = [H[1]/2, (H[1]+H[2]/2), (H[1]+H[2]+H[3])/2]
    mq = maximum(abs.(Qy))
    p = Plots.plot(Qy,-z,linewidth=2,ylabel="z [m]",xlabel=L"Q_y",ylimit=[-H[end],0],xlimit=[-mq,mq],
    xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(400,1000),margin=5mm)

    local savename = plotpath*"../Qy.png"
    savefig(p,savename)
end

function plot_unstable_vert(H,max_eve,psi_emp,plotpath,plotname,ell)
    max_eve_norm = max_eve./maximum(max_eve)
    z = [H[1]/2, (H[1]+H[2]/2), (H[1]+H[2]+H[3])/2]
    p = Plots.plot(max_eve_norm,-z,linewidth=2,ylabel="z [m]",xlabel=L"\psi_\mathrm{norm}",xlimit=[0,1.],ylimit=[-H[end],0],
    label=L"\psi_\mathrm{max} \ (L.S.)",xtickfontsize=12,ytickfontsize=12,xlabelfontsize=18,ylabelfontsize=18,legendfontsize=12,size=(400,1000),margin=5mm)
    Plots.plot!(psi_emp,-z,linewidth=2,label=L"\psi_\mathrm{max} \ (model)")

    local savename = @sprintf("%s_%04d.png", joinpath(plotpath, plotname), ell)
    savefig(p,savename)

end


