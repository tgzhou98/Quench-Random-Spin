using DelimitedFiles
using Statistics
using Plots
using LaTeXStrings
ENV["GKSwstype"]="100"
gr()

function convert_to_jld()
    integration_data = convert(Matrix, CSV.read("IntegrationData.csv"; header=false))
    abscissae1 = integration_data[1,:]
    weights1 = integration_data[2,:]
    errweights1 = integration_data[3,:]
    JLD.save("IntegrationData.jld", 
        "abscissae1", abscissae1,
        "weights1", weights1,
        "errweights1", errweights1,
        )
end



#################################################################
# hack Gadfly
#################################################################
function _findnz_val(x_vec::AbstractVector, y_vec::AbstractVector, testf::Function, A::AbstractMatrix{T}) where T
    N = Base.count(testf, A)
    is = zeros(Float64, N)
    js = zeros(Float64, N)
    zs = Array{T}(undef, N)
    if N == 0
        return (is, js, zs)
    end
    count = 1
    for j=1:size(A,2), i=1:size(A,1)
        Aij = A[i,j]
        if testf(Aij)
            is[count] = x_vec[j]
            js[count] = y_vec[i]
            zs[count] = Aij
            count += 1
        end
    end
    return (is, js, zs)
end


# function spy_x_y(x_vec::AbstractVector, y_vec::AbstractVector, M::AbstractMatrix, yflip_val::Bool = false, elements...; mapping...)
#     is, js, values = _findnz_val(x_vec, y_vec, x->!isnan(x), M)
#     # Assume that is and js is equal discreted
#     x_delta = x_vec[2] - x_vec[1]
#     y_delta = y_vec[2] - y_vec[1]
#     n,m = size(M)
#     Gadfly.plot(x=is, y=js, color=values,
#         Coord.cartesian(yflip=yflip_val, fixed=true, 
#         xmin=x_vec[1] - x_delta * 0.5, xmax=x_vec[m] + x_delta * 0.5,
#         ymin=y_vec[1] - y_delta * 0.5, ymax=y_vec[n] + y_delta * 0.5),
#         Scale.color_continuous,
#         Geom.rectbin,
#         Scale.x_continuous,
#         Scale.y_continuous,
#         elements...; mapping...)
# end

#################################################################
# Equilibrium part
#################################################################

function plot_ρω(ρω, abscissaeω, βJ, hinJ, field, Δ, β, plot_t::Float64 = 0.0, index::Int64 = 0, mode::String = "equilbrium")

    nω = length(abscissaeω)
    ρω_new = Array{Float64, 2}(undef, nω, 4)
    for a in 1:2
        for b in 1:2
            ρω_new[:, (a-1)*2+(b)] = real(ρω[:,a,b])
        end
    end
    beta_str = round(β, digits=5)
    if mode == "equilbrium"
        plot_t = 0.0
    end
    
    label_rhoomega = Array{String, 2}(undef, 1, 4)
    temp_index = 0
    for a in ["u", "d"]
        for b in ["u", "d"]
            temp_index += 1
            label_rhoomega[1, temp_index] = "\$ \\rho_{$(a) $(b)} \$"
        end
    end

    J = βJ / β
    
    myplot = Plots.plot(abscissaeω/J, ρω_new,
        title = latexstring("\$ \\textrm{t=}$(round(plot_t, digits=5)), \\beta=$(round(β, digits=5)), \\beta J=$(round(βJ, digits=5)), h/J=$(hinJ), f=$(Δ) \$"),
        label = label_rhoomega, xlabel=L"$ \omega/J $", ylabel=L"$ \rho (\omega) $", dpi=300)
    savefig(myplot, "ρω$(index)_$(mode).png")
    writedlm("ρω$(index)_$(mode).csv", hcat(abscissaeω/J, ρω_new), ',')
end


function plot_GRω(GRω, abscissaeω, βJ, hinJ, field, Δ, β, plot_t::Float64 = 0.0, index::Int64 = 0, mode::String = "equilbrium")

    nω = length(abscissaeω)
    GRω_new = Array{Float64, 2}(undef, nω, 4)
    a = 1
    for b in 1:2
        GRω_new[:, (a-1)*2+(b)] = real(GRω[:,1,b])
    end
    a = 2
    for b in 1:2
        GRω_new[:, (a-1)*2+(b)] = imag(GRω[:,1,b])
    end
    beta_str = round(β, digits=5)
    if mode == "equilbrium"
        plot_t = 0.0
    end
    
    label_GRomega = Array{String, 2}(undef, 1, 4)
    temp_index = 0
    for a in ["\\textrm{Re}", "\\textrm{Im}"]
        for b in ["u u", "u d"]
            temp_index += 1
            label_GRomega[1, temp_index] = "\$ $a  G_{$(b)}^R(\\omega) \$"
        end
    end

    J = βJ / β
    
    myplot = Plots.plot(abscissaeω/J, GRω_new,
        title = latexstring("\$ \\textrm{t=}$(round(plot_t, digits=5)), \\beta=$(round(β, digits=5)), \\beta J=$(round(βJ, digits=5)), h/J=$(hinJ), f=$(Δ) \$"),
        label = label_GRomega, xlabel=L"$ \omega/J $", ylabel=L"$ G^R (\omega) $", dpi=300)
    savefig(myplot, "GRω$(index)_$(mode).png")
    writedlm("GRω$(index)_$(mode).csv", hcat(abscissaeω/J, GRω_new), ',')
end

function save_arrayω(filename, arrayomega, abscissaeω, real_or_imag, index::Int64 = 0, mode::String = "equilbrium")

    nω = length(abscissaeω)
    arrayomega_new = Array{Float64, 2}(undef, nω, 4)

    for a in 1:2
        for b in 1:2
            if real_or_imag == "real"
                arrayomega_new[:, (a-1)*2+(b)] = real(arrayomega[:,a,b])
            elseif real_or_imag == "imag"
                arrayomega_new[:, (a-1)*2+(b)] = imag(arrayomega[:,a,b])
            end
        end
    end

    writedlm("$(filename)_$(real_or_imag)$(index)_$(mode).csv", arrayomega_new, ',')
end


function plot_Gt_eq(Gt, abscissaet, βJ, hinJ, field, Δ, β, plot_t::Float64 = 0.0, index::Int64 = 0, mode::String = "equilbrium")

    nt = length(abscissaet)
    Gt_new = Array{Float64, 2}(undef, nt, 4)
    a = 1
    for b in 1:2
        Gt_new[:, (a-1)*2+(b)] = real(Gt[:,1,b])
    end
    a = 2
    for b in 1:2
        Gt_new[:, (a-1)*2+(b)] = imag(Gt[:,1,b])
    end

    # beta_str = round(β, digits=5)
    if mode == "equilbrium"
        plot_t = 0.0
    end

    label_Gt = Array{String, 2}(undef, 1, 4)
    temp_index = 0
    
    for a in ["\\textrm{Re}", "\\textrm{Im}"]
        for b in ["u u", "u d"]
            temp_index += 1
            label_Gt[1, temp_index] = "\$ $a  G_{$(b)}(t) \$"
        end
    end
    
    myplot = Plots.plot(abscissaet/β, Gt_new,
        title = latexstring("\$ \\textrm{t=}$(round(plot_t, digits=5)), \\beta=$(round(β, digits=5)), \\beta J=$(round(βJ, digits=5)), h/J=$(hinJ), f=$(Δ) \$"),
        label = label_Gt, xlabel = L"$ t/ \beta $", ylabel = L"$ G_<(t) $", dpi=300)
    savefig(myplot, "Gt$(index)_$(mode).png")
    writedlm("Gt$(index)_$(mode).csv", hcat(abscissaet/β, Gt_new), ',')
    
end



#################################################################
# phase diagram part
#################################################################
function plot_phase_diagram(betaJlist::Array{Float64, 1}, hlist::Array{Float64, 1}, ρω_centerlist::Array{Float64, 2},
    Sxlist::Array{Float64,2}, h0list::Array{Float64,2}, fix_beta::Bool)

    h_num = length(hlist)
    # hlist = hlist / 5.0
    if fix_beta
        xstr = "\$h \\beta\$"
        ystr = "\$\\beta J, \\ \\beta =2\\pi \$"
        hlist = round.(hlist*2π, digits=2)
    else
        xstr = "\$h / J \$"
        ystr = "\$\\beta J, \\ J =5/(2\\pi) \$"
    end

    aspect_ratio_val = 1/20
    

    myplot = Plots.contour(hlist, betaJlist, ρω_centerlist, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } \\rho(\\omega=0)\$", dpi=300)
    savefig(myplot, "phase_diagram_rhoomega_contour.png")

    myplot = Plots.contour(hlist, betaJlist, Sxlist, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } S_x\$", dpi=300)
    savefig(myplot,"phase_diagram_Sx.png")

    myplot = Plots.contour(hlist, betaJlist, h0list, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } h_0\$", dpi=300)
    savefig(myplot,"phase_diagram_h0_contour.png")


    myplot = Plots.heatmap(hlist, betaJlist, ρω_centerlist, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } \\rho(\\omega=0)\$", dpi=300)
    savefig(myplot, "phase_diagram_rhoomega_heatmap.png")

    myplot = Plots.heatmap(hlist, betaJlist, Sxlist, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } S_x\$", dpi=300)
    savefig(myplot,"phase_diagram_Sx_heatmap.png")

    myplot = Plots.heatmap(hlist, betaJlist, h0list, aspect_ratio=aspect_ratio_val,
        xlims = (hlist[1],hlist[h_num]), xlabel=xstr, ylabel=ystr, title = "\$\\textrm{Phase diagram of } h_0\$")
    savefig(myplot,"phase_diagram_h0_heatmap.png")

    #########################################################################################

    # myplot = Plots.contour(hlist, betaJlist, ρω_centerlist, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of center point of } \\rho\\omega\$")
    # savefig(myplot, "phase_diagram_rhoomega_contour.pdf")

    # myplot = Plots.contour(hlist, betaJlist, Sxlist, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of } S_x\$")
    # savefig(myplot,"phase_diagram_Sx.pdf")

    # myplot = Plots.contour(hlist, betaJlist, h0list, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of } h_0\$")
    # savefig(myplot,"phase_diagram_h0_contour.pdf")


    # myplot = Plots.heatmap(hlist, betaJlist, ρω_centerlist, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of center point of } \\rho\\omega\$")
    # savefig(myplot, "phase_diagram_rhoomega_heatmap.pdf")

    # myplot = Plots.heatmap(hlist, betaJlist, Sxlist, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of } S_x\$")
    # savefig(myplot,"phase_diagram_Sx_heatmap.pdf")

    # myplot = Plots.heatmap(hlist, betaJlist, h0list, aspect_ratio=1,
    #     xlims = (hlist[1],hlist[h_num]), xlabel="\$h J\$", ylabel=ystr, title = "\$\\textrm{Phase diagram of } h_0\$")
    # savefig(myplot,"phase_diagram_h0_heatmap.pdf")

end

#################################################################
# SY_dynamics part
#################################################################

function plot_Gt(GDynamicsttp, βJ, hinJ, field, Δ, β, begin_val, end_val, Δt, center_pt)

    # myplot = spy(real(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,1]), Guide.yticks(ticks=(collect(begin_val:end_val).-center_pt)*Δt),
    # Guide.xticks(ticks=(collect(begin_val:end_val).-center_pt)*Δt))
    # img = PNG("GLessertreal_upup.png", 6inch, 4inch,dpi=300)
    # draw(img, myplot)
    time_list = (collect(begin_val:end_val).-center_pt)*Δt / β
    

    xstr = "\$ t_1 / \\beta \$"
    ystr = "\$ t_2 / \\beta \$"


    # myplot = Plots.contour(time_list, time_list, real(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,1]), aspect_ratio=1,
    #     xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Re} G^<_{11}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    # savefig(myplot, "GLessertreal_upup_contour.pdf")

    myplot = Plots.heatmap(time_list, time_list, real(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,1]), aspect_ratio=1, dpi=600,
        xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Re} G^<_{11}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    savefig(myplot, "GLessertreal_upup_heapmap.png")

    # myplot = Plots.contour(time_list, time_list, imag(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,1]), aspect_ratio=1,
    #     xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Im} G^<_{11}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    # savefig(myplot, "GLessertimag_upup_contour.pdf")

    myplot = Plots.heatmap(time_list, time_list, imag(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,1]), aspect_ratio=1, dpi=600,
        xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Im} G^<_{11}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    savefig(myplot, "GLessertimag_upup_heapmap.png")



    # myplot = Plots.contour(time_list, time_list, real(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,2]), aspect_ratio=1,
    #     xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Re} G^<_{12}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    # savefig(myplot, "GLessertreal_updown_contour.pdf")

    myplot = Plots.heatmap(time_list, time_list, real(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,2]), aspect_ratio=1, dpi=600,
        xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Re} G^<_{12}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    savefig(myplot, "GLessertreal_updown_heapmap.png")

    # myplot = Plots.contour(time_list, time_list, imag(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,2]), aspect_ratio=1,
    #     xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Im} G^<_{12}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    # savefig(myplot, "GLessertimag_updown_contour.pdf")

    myplot = Plots.heatmap(time_list, time_list, imag(GDynamicsttp[begin_val:end_val,begin_val:end_val,1,2]), aspect_ratio=1, dpi=600,
        xlims = (time_list[1],time_list[length(time_list)]), xlabel=xstr, ylabel=ystr, title = "\$\\mathrm{Im} G^<_{12}(t_1, t_2) \\textrm{ Evolvolution ,} \\Delta t = $(Δt) \$")
    savefig(myplot, "GLessertimag_updown_heapmap.png")


    # begin_val = 990
    # # end_val = 1100
    # end_val = 1030
    # df = DataFrame()
    # for a in 1:2
    #     for b in 1:2
    #         df=vcat(df, DataFrame(x=(collect(begin_val:end_val).-1001)*50/1000, y=imag(diag(GGreaterDynamicsttp[:,:,a,b], -1)[begin_val:end_val]), label="ρ$(a)$(b)"))
    #     end
    # end

    # myplot = Gadfly.Gadfly.plot(df, x="x", y="y", color="label", Geom.line, 
    #     Guide.xlabel("Time"), Guide.ylabel("Im(G(t, t-1))"), Guide.title("G(t, t-1) βJ$(βJ), h$(hinJ)J, $(field), Delta$(Δ)"),
    #     Scale.discrete_color_manual("blue","red", "green"))
    # img = PDF("GGreatertimag.pdf", 6inch, 4inch)
    # draw(img, myplot)
end


function plot_Sx(GDynamicsttp, βJ, hinJ, field, Δ, β, begin_val, end_val, Δt, center_pt)
    # Assume a pass the GLesserDynamicsttp
    # S_x = Im( i (<c_{\uparrow}^{\dagger}(t') c_{\downarrow}(t)> + <c_{\downarrow}^{\dagger}(t') c_{\uparrow}(t)>)/2 )
    #     = Im (G<_{\uparrow, \downarrow}(t, t'))

    beta_str = round(β, digits=5)
    a=1
    b=2
    
    myplot = Plots.plot((collect(begin_val:end_val).-center_pt)*Δt / β,
        imag(diag(GDynamicsttp[:,:,a,b])[begin_val:end_val]),
        title = latexstring("\$ \\beta=$(beta_str), \\beta J=$(βJ), h/J=$(hinJ), f=$(Δ) \$"),
        xlabel = L"$ t / \beta $", ylabel = L"$ S_x(t) $", dpi=300)
    savefig(myplot, "Sx_updown.png")
end


function plot_allS(GDynamicsttp, βJ, hinJ, JbarinJ, field, Δ, β, begin_val, end_val, Δt, center_pt)
    # Assume a pass the GLesserDynamicsttp
    # S_x = Im( i (<c_{\uparrow}^{\dagger}(t') c_{\downarrow}(t)> + <c_{\downarrow}^{\dagger}(t') c_{\uparrow}(t)>)/2 )
    #     = Im (G<_{\uparrow, \downarrow}(t, t'))

    beta_str = round(β, digits=5)
    # a=1
    # b=2


    Sall_new = Array{Float64, 2}(undef, end_val - begin_val + 1, 4)

    Sall_new[:, 1] = real(1/2.0*(-1.0im)*(diag(GDynamicsttp[:,:,1,2])[begin_val:end_val] + diag(GDynamicsttp[:,:,2,1])[begin_val:end_val]))
    Sall_new[:, 2] = real(1/2.0*(diag(GDynamicsttp[:,:,1,2])[begin_val:end_val] - diag(GDynamicsttp[:,:,2,1])[begin_val:end_val]))
    Sall_new[:, 3] = real(1/2.0*(-1.0im)*(diag(GDynamicsttp[:,:,1,1])[begin_val:end_val] - diag(GDynamicsttp[:,:,2,2])[begin_val:end_val]))
    Sall_new[:, 4] = sqrt.(Sall_new[:, 1].^2 + Sall_new[:, 2].^2 + Sall_new[:, 3].^2) 


    label_St = Array{String, 2}(undef, 1, 4)
    label_St[1] = L"$ S_x(t) $"
    label_St[2] = L"$ S_y(t) $"
    label_St[3] = L"$ S_z(t) $"
    label_St[4] = L"$ \sqrt(S^2)(t) $"
    
    myplot = Plots.plot((collect(begin_val:end_val).-center_pt)*Δt / β,
        Sall_new,
        title = latexstring("\$ \\bar{J}/J =$(JbarinJ), \\beta=$(beta_str), \\beta J=$(βJ), h/J=$(hinJ), f=$(Δ) \$"),
        label = label_St,
        xlabel = L"$ t / \beta $", ylabel = L"$ S_i(t) $", dpi=300)
    savefig(myplot, "S_all.png")

    writedlm("S_all.csv", hcat((collect(begin_val:end_val).-center_pt)*Δt / β, Sall_new), ',')
end


function plot_n_up(GDynamicsttp, βJ, hinJ, field, Δ, β, begin_val, end_val, Δt, center_pt)
    # Assume a pass the GLesserDynamicsttp
    # n_up = Im( i (<c_{\uparrow}^{\dagger}(t') c_{\uparrow}(t)>) )
    #     = Im (G<_{\uparrow, \uparrow}(t, t'))
    beta_str = round(β, digits=5)
    a=1
    b=1
    
    myplot = Plots.plot((collect(begin_val:end_val).-center_pt)*Δt / β,
        imag(diag(GDynamicsttp[:,:,a,b])[begin_val:end_val]),
        title = latexstring("\$ \\beta=$(beta_str), \\beta J=$(βJ), h/J=$(hinJ), f=$(Δ) \$"),
        xlabel = L"$ t / \beta $", ylabel = L"$ n_{\uparrow \uparrow}(t) $", dpi=300)
    savefig(myplot, "n_up.png")
end

function export_Sx(GDynamicsttp)
    # Assume a pass the GLesserDynamicsttp
    # i <c^{\dagger}(t') c(t)>
    # writedlm("Sx_$(βJ)_$(hinJ)_$(field)_$(Δ).csv",  imag(diag(GDynamicsttp[:,:,1,2])), ',')
    writedlm("Sx.csv",  imag(diag(GDynamicsttp[:,:,1,2])), ',')
end


function export_cutoff(path_name::String, Λω::Float64, Λt_arr::Array{Float64, 2}, time::Int64)
    # Assume a pass the GLesserDynamicsttp
    # i <c^{\dagger}(t') c(t)>
    # writedlm("Sx_$(βJ)_$(hinJ)_$(field)_$(Δ).csv",  imag(diag(GDynamicsttp[:,:,1,2])), ',')
    writedlm("$(path_name)/deltaT.csv",  Λt_arr / time, ',')
    writedlm("$(path_name)/CutoffT.csv",  Λt_arr, ',')
    writedlm("$(path_name)/Cutoffomega.csv",  Λω, ',')
end


function findlocalmaxima(signal::Vector)
    # Some ad hoc algorithm
    inds = Int[]
    max_val = signal[1]
    max_pos = 1
    flag = 1
    if length(signal)>1
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]
                max_val = signal[i]
                max_pos = i
                flag = 1
            end
            if flag == 1 && signal[max_pos] > 0.5 && (-signal[i] + signal[max_pos])/signal[max_pos] > 0.5
                push!(inds, max_pos)
                flag = 2
            end
            # Assume the base line is 0
        end
    end
    inds
  end

function peak_analyze(run_num::Int64, ρωlist::Array{ComplexF64, 4}, abscissaeω::Array{Float64, 1}, betaJlist::Array{Float64, 1})
    for i in 1:run_num
        ρω = ρωlist[i,:,:,:]
        peaks_pt = findlocalmaxima(real(ρω[:,1,1]))
        peaks_ω = [abscissaeω[i] for i in peaks_pt] 
        # println(peaks)
        # println([real(ρω[i,1,1]) for i in peaks])
        @assert length(peaks_pt)==4
        nω = size(ρω)[1]
        center_pt = (nω+1)÷2
        ralative_peaks = abs.(peaks_ω .- abscissaeω[center_pt])
        ratio_peaks12 = (ralative_peaks[1]/ralative_peaks[2] + ralative_peaks[4]/ralative_peaks[3])/2
        println("βJ=$(betaJlist[i]), peaks ratio = $(ratio_peaks12)")
    end
end

function spectral_ratio(ρω::Array{ComplexF64, 3})
    rho_upup = real.(ρω[:,1,1])

    nω = length(rho_upup)
    rho_peak = maximum(rho_upup)
    rho_0 = rho_upup[(nω+1)÷2]
    
    ratio_val = rho_0 / rho_peak

    return ratio_val
end

