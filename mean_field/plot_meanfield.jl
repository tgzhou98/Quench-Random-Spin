using Plots
using LinearAlgebra
using LaTeXStrings
using Distributions
using JLD2
using LsqFit
using DelimitedFiles




function A_polynomial(xi::Array{Float64, 1})
    A = (-xi[1]^2+xi[2]^2-4*xi[2]*xi[3]+xi[3]^2)
    omega = 0.0
    if A > 0
        omega = sqrt(A)
    else
        omega = 0.0
    end

    return omega
end

function B_polynomial(xi::Array{Float64, 1})
    B = (xi[1]^2+xi[2]^2+xi[3]^2)
    decay = 0.0

    A = (-xi[1]^2+xi[2]^2-4*xi[2]*xi[3]+xi[3]^2)
    if A > 0
        decay = sqrt(B)
    else
        decay = -sqrt(-A) + sqrt(B)
    end

    return decay
end

# xi_arr_all = [[1.0, 0.8, -1.5], [1.0, 0.8, 1.5], [1.0, 1.0, -2.0], [1.0, 1.0, 2.0]]

# xi_arr_all = [[-1.46, 1.89, -1.25], [1.47, -1.78, 
#     0.27], [0.36, -1.66, 0.88], [1.79, -1.68, 0.99], [-0.8, 
#     1.07, -0.76], [-0.51, 0.97, -1.51], [0.25, 
#     0.08, -1.42], [-1.58, 1.05, -1.88], [-0.01, -0.29, 
#     1.34], [0.09, -1.54, 0.18], [1.51, -0.22, 1.54], [0.25, 
#     0.06, 0.81], [-1.49, -1.27, 
#     0.83], [0.89, -0.1, -1.66], [-0.73, 1.89, 0.22], [0.56, 
#     0.58, 0.93], [-1.83, -0.97, -1.04], [-1.15, 1.23, 
#     1.27], [1.31, -1.68, -1.31], [1.2, 1.07, 0.]]

# xi_arr_all = []
# xi_xarr = collect(-0.4:0.05:0.2)
# for xi_x in xi_xarr
#     push!(xi_arr_all, [0.2,xi_x,-0.2-xi_x])
# end
# xi_arr_all = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 0.2], [-0.05, -0.15, 0.2], [0.0, -0.2, 0.2], [0.05, -0.25, 0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]]

J = 1.0
# tspan = (0.0,12.0)
tspan = (0.0,16.0)
delta_t = 0.04
t_arr_fix = collect(tspan[1]:(delta_t):tspan[2])
t_num_set = length(t_arr_fix)
nsite = 2000
n_alpha = 16
beta_h = 0.4
Nrealization = 20
fontsize = 11

# rhox
# xi_arr_all = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 0.2], [-0.05, -0.15, 0.2], [0.0, -0.2, 0.2], [0.05, -0.25, 0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]]
# rhoy
# xi_arr_all = [[0.2, 0.2, -0.4], [0.15, 0.2, -0.35], [0.1, 0.2, -0.3], [0.05, 0.2, -0.25], [0.0, 0.2, -0.2], [-0.05, 0.2, -0.15], [-0.1, 0.2, -0.1], [-0.15, 0.2, -0.05], [-0.2, 0.2, 0.0], [-0.25, 0.2, 0.05], [-0.3, 0.2, 0.1], [-0.35, 0.2, 0.15], [-0.4, 0.2, 0.2]]
# rhoz
# xi_arr_all = [[0.2, -0.4, 0.2], [0.2, -0.35, 0.15], [0.2, -0.3, 0.1], [0.2, -0.25, 0.05], [0.2, -0.2, 0.0], [0.2, -0.15, -0.05], [0.2, -0.1, -0.1], [0.2, -0.05, -0.15], [0.2, 0.0, -0.2], [0.2, 0.05, -0.25], [0.2, 0.1, -0.3], [0.2, 0.15, -0.35], [0.2, 0.2, -0.4]]

direct_string = ["x", "y", "z"]


xi_arr_all_x = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 0.2], [-0.05, -0.15, 0.2], [0.0, -0.2, 0.2], [0.05, -0.25, 0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]]
xi_arr_all_y = [[0.2, 0.2, -0.4], [0.15, 0.2, -0.35], [0.1, 0.2, -0.3], [0.05, 0.2, -0.25], [0.0, 0.2, -0.2], [-0.05, 0.2, -0.15], [-0.1, 0.2, -0.1], [-0.15, 0.2, -0.05], [-0.2, 0.2, 0.0], [-0.25, 0.2, 0.05], [-0.3, 0.2, 0.1], [-0.35, 0.2, 0.15], [-0.4, 0.2, 0.2]]
xi_arr_all_z = [[0.2, -0.4, 0.2], [0.2, -0.35, 0.15], [0.2, -0.3, 0.1], [0.2, -0.25, 0.05], [0.2, -0.2, 0.0], [0.2, -0.15, -0.05], [0.2, -0.1, -0.1], [0.2, -0.05, -0.15], [0.2, 0.0, -0.2], [0.2, 0.05, -0.25], [0.2, 0.1, -0.3], [0.2, 0.15, -0.35], [0.2, 0.2, -0.4]]


xi_x_len = length(xi_arr_all_x)

omega_fit_mean_error_xyz = zeros(xi_x_len, 3, length(direct_string))
decay_fit_mean_error_xyz = zeros(xi_x_len, 3, length(direct_string))
for adirect in 1:3
    # adirect = 2

    xi_arr_all = xi_arr_all_x
    if adirect == 1
        xi_arr_all = xi_arr_all_x
    elseif adirect == 2
        xi_arr_all = xi_arr_all_y
    elseif adirect == 3
        xi_arr_all = xi_arr_all_z
    end

    # scaling_c = [0.951, 0.948, 0.899]
    scaling_c = [0.9876, 0.9869, 0.9054]

    xi_len = length(xi_arr_all)
    Sxyz_t_aver_xi_all = zeros(xi_len, t_num_set, 3)
    total_dirname = "figs_exp_$(direct_string[adirect])"
    # fitstep_span = (1, 12*25)
    # fitstep_span = (2, 12*25)

    if !isdir(total_dirname)
        mkdir(total_dirname)
    end
    for (i,xi_arr) in enumerate(xi_arr_all)
        #  tJ_arr_fix, Sxyz_t_aver = main(nsite, n_alpha, beta_h, Nrealization, J, xi_arr, tspan, delta_t)

        fileprefix = "xi$(xi_arr)_beta_h$(beta_h)_Nrand$(Nrealization)_tspan$(tspan)"
        # total_dirname = "figs"
        Sxyz_t_aver = JLD2.load("$(total_dirname)/Sxyz_meanfield_$(fileprefix).jld2", "Sxyz_t_aver")
        Sxyz_t_aver_xi_all[i,:,:] = Sxyz_t_aver
    end


    fit_end_arr = collect(320:10:400)
    omega_diff_fit = zeros(xi_x_len, length(fit_end_arr))
    decay_diff_fit = zeros(xi_x_len, length(fit_end_arr))
    for (iind, fit_end) in enumerate(fit_end_arr)
        fitstep_span = (10, fit_end)

        coef_fit_arr = zeros(xi_x_len, 5)
        for ll in 1:xi_x_len
            println("ll=$(ll), iind=$(fit_end)")
            # @. model(x, p) = p[1]*exp(-x*p[2])*cos(p[3]*x+p[4]) + p[5]
            @. model(x, p) = p[1]*1.0/cosh(x*p[2])*cos(p[3]*x+p[4]) + p[5]
            xdata = (t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]
            ydata = (Sxyz_t_aver_xi_all[ll,:,1])[fitstep_span[1]:fitstep_span[2]]
            p0 = [0.516627,  2.240574,    0.31,    -0.728552, 0.0]
            lb = [0.2,0.0,0.0,-π, -0.01]
            ub = [2.0,4.0,2.0,π, 0.01]
            fitfunc = curve_fit(model, xdata, ydata, p0, lower=lb, upper=ub)
            # fitfunc = curve_fit(model, xdata, ydata, p0)

            coef_fit_arr[ll,:] = coef(fitfunc)
        end
        omega_diff_fit[:,iind] = coef_fit_arr[:,3]
        decay_diff_fit[:,iind] = coef_fit_arr[:,2]
    end

    println("omega_diff_fit")
    display(omega_diff_fit)
    for ll in 1:xi_x_len
        omega_fit_mean_error_xyz[ll,1,adirect] = mean(omega_diff_fit[ll,:])
        omega_fit_mean_error_xyz[ll,2,adirect] = minimum(omega_diff_fit[ll,:]) - omega_fit_mean_error_xyz[ll,1,adirect]
        omega_fit_mean_error_xyz[ll,3,adirect] = maximum(omega_diff_fit[ll,:]) - omega_fit_mean_error_xyz[ll,1,adirect]

        decay_fit_mean_error_xyz[ll,1,adirect] = mean(decay_diff_fit[ll,:])
        decay_fit_mean_error_xyz[ll,2,adirect] = minimum(decay_diff_fit[ll,:]) - decay_fit_mean_error_xyz[ll,1,adirect]
        decay_fit_mean_error_xyz[ll,3,adirect] = maximum(decay_diff_fit[ll,:]) - decay_fit_mean_error_xyz[ll,1,adirect]
    end


    # a two-parameter exponential model
    # x: array of independent variables
    # p: array of model parameters
    # model(x, p) will accept the full data set as the first argument `x`.
    # This means that we need to write our model function so it applies
    # the model to the full dataset. We use `@.` to apply the calculations
    # across all rows.
    coef_fit_arr = zeros(length(xi_arr_all), 5)
    # fitstep_span = (3, 120)
    fitstep_span = (1, 16*25)
    for ll in 1:length(xi_arr_all)
        # @. model(x, p) = p[1]*exp(-x*p[2])*cos(p[3]*x+p[4]) + p[5]
        @. model(x, p) = p[1]*1.0/cosh(x*p[2])*cos(p[3]*x+p[4]) + p[5]
        xdata = (t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]
        ydata = (Sxyz_t_aver_xi_all[ll,:,1])[fitstep_span[1]:fitstep_span[2]]
        p0 = [0.516627,  0.540574,    0.71,    -0.728552, 0.0]
        lb = [0.2,0.0,0.0,-π, -0.01]
        ub = [2.0,4.0,2.0,π, 0.01]
        fitfunc = curve_fit(model, xdata, ydata, p0, lower=lb, upper=ub)

        coef_fit_arr[ll,:] = coef(fitfunc)
        display(coef_fit_arr)
    end

    fit_all_xi = zeros(length(xi_arr_all), (fitstep_span[2]-fitstep_span[1]+1))
    for ll in 1:length(xi_arr_all)
        # @. modelfit(x, p) = p[1]*exp(-x*p[2])*cos(p[3]*x+p[4]) + p[5]
        @. modelfit(x, p) = p[1]*1.0/cosh(x*p[2])*cos(p[3]*x+p[4]) + p[5]
        pval = coef_fit_arr[ll,:]
        xdata = (t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]
        for (i,x) in enumerate(xdata)
            a = modelfit(x, pval)
            fit_all_xi[ll,i] = a
        end
    end


    label_xi = Array{String, 2}(undef, 1, length(xi_arr_all))
    for ll in 1:length(xi_arr_all)
        label_xi[1,ll] = "\$ \\xi=$(xi_arr_all[ll]) \$"
    end
    plot(t_arr_fix*J, transpose(Sxyz_t_aver_xi_all[:,:,1]), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization), \rho_{%$(direct_string[adirect])}$", xlabel=L"$tJ$",ylabel=L"\langle\hat{S}^{x}(t)\rangle", legend = false, titlefontsize=fontsize, legendfont=font(fontsize))
    figplot = plot!(((t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]), transpose(fit_all_xi), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization), \rho_{%$(direct_string[adirect])}$", xlabel=L"$tJ$",ylabel=L"\langle\hat{S}^{x}(t)\rangle", legend = false, titlefontsize=fontsize, legendfont=font(fontsize), linestyle=:dash)

    fileprefix = "beta_h$(beta_h)_Nrand$(Nrealization)_tspan$(tspan)"

    println("$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect])")
    savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).pdf")
    savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).png")

    println(coef_fit_arr[:,3])

    # Sxfreqtheory = zeros(length(xi_arr_all))
    Sxfreqtheory = [A_polynomial(xi) for xi in xi_arr_all]
    Sxdecaytheory = [B_polynomial(xi) for xi in xi_arr_all]
    # SxfreqKB = [3.34044, 1.89207, 2.72674, 2.75744, 1.9145, 2.66232, 1.4702, 3.09645, 
    # 1.68609, 1.75667, 1.42522, 0.602516, 2.10143, 1.27365, 1.36014, 
    # 0.0070625, 0.000278918, 0.000200615, 0.000639637, 0.350797]

    xix_arr = collect(-0.4:0.05:0.2)
    # scatter(xix_arr, hcat(coef_fit_arr[:,3], Sxfreqtheory), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$\mathrm{Cases}$",ylabel=L"\Omega_x", label=["meanfield" "theory" "Large N"],  titlefontsize=fontsize, legendfont=font(fontsize))

    # @. model_freq(x, p) = p[1]*real(sqrt(complex(x+p[2])))
    # xdata = Sxfreqtheory.^2
    # ydata = coef_fit_arr[:,3]
    # p0 = [1.0, 0.1]
    # fitsqrtAfunc = curve_fit(model_freq, xdata, ydata, p0)

    # p0val = coef(fitsqrtAfunc)

    # fitfreq = zeros(length(xdata))
    # for (i,x) in enumerate(xdata)
    #     a = model_freq(x, p0val)
    #     fitfreq[i] = a
    # end

    # omega_list_xyz = [0.0 0.6263006677264809 0.5740071344287895; 0.0 0.5889913577181636 0.41381231158131015; 0.0 0.5394047662598402 0.2450183703916318; 5.6649544275805435e-9 0.4862911382984757 0.07527734173711145; 4.2873988391116597e-8 0.43284142215869875 0.016796370357462787; 0.07343244215565856 0.3796713740167248 0.012061032207967082; 0.31296244340860657 0.31173313296665073 0.0; 0.3792967466313909 0.09231088301946341 2.7027748707273832e-8; 0.4314399212449717 6.594175089257451e-8 0.0; 0.4859903621423288 1.5315567249590236e-8 0.06700920500306148; 0.5419444898790315 0.0 0.25547310051140654; 0.5944236155711664 0.0 0.4356040960487336; 0.6333699799987746 0.0 0.5951850043610303]
    omega_list_xyz = [1.4505593718982169e-9 0.6298179888873356 0.5800281711970151; 1.60478596014347e-8 0.588664650024305 0.4362877024785485; 4.8519359584529706e-8 0.5383628302200096 0.2868370994097614; 2.2400788050567584e-8 0.48541541679735195 0.16520049473146733; 0.03143998115292887 0.4322156184943408 0.05905719995372856; 0.16717408480406637 0.3789100461125505 0.0015230681110718724; 0.31636478738125734 0.3142897179850306 3.715195092723649e-8; 0.37893321104799976 0.17007417703616665 0.0; 0.43087304638143725 0.04804582737819918 0.04936731250692237; 0.48515935438442104 2.6268597923794813e-7 0.16616638252795274; 0.540777851603687 1.3767544330792227e-7 0.2969490073247213; 0.5938729105580676 1.4651421632314675e-8 0.45740195501487774; 0.6371055671561815 0.0 0.5999067413904184]

    # scatter(xix_arr, hcat(coef_fit_arr[:,3], Sxfreqtheory, fitfreq), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$\xi_x$",ylabel=L"\Omega_x", label=["mean field" "\$J\\sqrt{A}\$" "\$$(p0val[1])J\\sqrt{A$(p0val[2])}\$"],  titlefontsize=fontsize, legendfont=font(fontsize))
    dataplot = hcat(coef_fit_arr[:,3], Sxfreqtheory*scaling_c[adirect], omega_list_xyz[:,adirect])
    scatter(xix_arr, dataplot, title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization), \rho_{%$(direct_string[adirect])}$", xlabel=L"$\xi_x$",ylabel=L"\omega/J", label=["New Theory" "\$ $(scaling_c[adirect]) J\\sqrt{A}\$" "Exp"],  titlefontsize=fontsize, legendfont=font(fontsize))
    savefig("$(total_dirname)/Sxfreq_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).pdf")
    savefig("$(total_dirname)/Sxfreq_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).png")

    savefig("fig_expdata/Sxfreq_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).pdf")
    savefig("fig_expdata/Sxfreq_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).png")


    datadecayplot = hcat(coef_fit_arr[:,2], Sxdecaytheory)
    scatter(xix_arr, datadecayplot, title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization), \rho_{%$(direct_string[adirect])}$", xlabel=L"$\xi_x$",ylabel=L"\Gamma/J", label=["New Theory" "\$ J\\sqrt{\\mathbf{\\xi}}\$" "Exp"],  titlefontsize=fontsize, legendfont=font(fontsize))
    savefig("$(total_dirname)/Sxdecay_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).pdf")
    savefig("$(total_dirname)/Sxdecay_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).png")
    savefig("fig_expdata/Sxdecay_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).pdf")
    savefig("fig_expdata/Sxdecay_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).png")


    JLD2.save("$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).jld2", "tJ_arr", t_arr_fix*J, "Sxyz_t_aver_xi_all", Sxyz_t_aver_xi_all, "coef_fit_arr", coef_fit_arr, "plot_data", dataplot)

    writedlm("fig_expdata/Sxyz_meanfield_$(fileprefix)_xiall_rho$(direct_string[adirect]).csv",  dataplot, ',')
    writedlm("fig_expdata/Sxyz_meanfield_decay_$(fileprefix)_xiall_rho$(direct_string[adirect]).csv",  datadecayplot, ',')
    
end


total_dirname = "meanfield_data"
# fitstep_span = (1, 12*25)

if !isdir(total_dirname)
     mkdir(total_dirname)
end
for a in 1:3
    writedlm("$(total_dirname)/meanfield_omega_meanerrorbar_rho$(direct_string[a]).csv",  omega_fit_mean_error_xyz[:,:,a], ',')
    writedlm("$(total_dirname)/meanfield_decay_meanerrorbar_rho$(direct_string[a]).csv",  decay_fit_mean_error_xyz[:,:,a], ',')
end