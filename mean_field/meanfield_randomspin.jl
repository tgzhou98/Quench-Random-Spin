using DifferentialEquations
using Plots
using LinearAlgebra
using LaTeXStrings
using Distributions
using JLD2

#=
f(u,p,t) = 1.01*(u)
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")
=#

# function parameterized_lorenz!(du,u,p,t)
#  du[1] = p[1]*(u[2]-u[1])
#  du[2] = u[1]*(p[2]-u[3]) - u[2]
#  du[3] = u[1]*u[2] - p[3]*u[3]
# end

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


function parameterized_randomspin!(du,u,p,t)
     plength = length(p)
     nsite = convert(Int64, p[plength])
     # Jij_num = nsite*(nsite+1)÷2
     Jij_num = nsite*(nsite)
     xi_arr = p[1:3]
     Jij_arr = p[4:Jij_num+3]

     # initial
     for j in 1:nsite
          for l in 1:3
               du[3*(j-1)+l] = 0.0
          end
     end

     for j in 1:nsite
          for i in 1:nsite
               du[3*(j-1)+1] += Jij_arr[(j-1)*nsite+i] * (xi_arr[2] * u[3*(i-1)+2] * u[3*(j-1)+3] - xi_arr[3] * u[3*(i-1)+3] * u[3*(j-1)+2])
               du[3*(j-1)+2] += Jij_arr[(j-1)*nsite+i] * (xi_arr[3] * u[3*(i-1)+3] * u[3*(j-1)+1] - xi_arr[1] * u[3*(i-1)+1] * u[3*(j-1)+3])
               du[3*(j-1)+3] += Jij_arr[(j-1)*nsite+i] * (xi_arr[1] * u[3*(i-1)+1] * u[3*(j-1)+2] - xi_arr[2] * u[3*(i-1)+2] * u[3*(j-1)+1])
          end
     end
end
function main(nsite, n_alpha, beta_h, Nrealization, J, xi_arr, tspan, delta_t, total_dirname)

     # J = 1.0
     # xi_arr = [1.0, 1.0, -0.2]
     # xi_arr = [1.0, 1.0, -0.5]
     # xi_arr = [1.0, 1.0, 2.0]
     # tspan = (0.0,50.0)

     # tspan = (0.0,12.0)
     # delta_t = 0.04
     # Jij_num = nsite*(nsite+1)÷2
     Jij_num = nsite*(nsite)
     sigmasqrt = sqrt(4*J^2/nsite/n_alpha)

     println("xi_arr=$(xi_arr), J=$(J), N=$(nsite), n_alpha=$(n_alpha)")


     # cubic spline interpolation
     # t_num_set = 401
     t_arr_fix = collect(tspan[1]:(delta_t):tspan[2])
     t_num_set = length(t_arr_fix)
     Sxyz_t_rand = zeros(t_num_set, 3, Nrealization)

     for irand in 1:Nrealization
          println("n_rand=$(irand), t_num=$(t_num_set)")

          # sx_arr = rand(Uniform(-W_disorder, W_disorder), nsite) .+ sinh(beta_h)
          # sy_arr = rand(Uniform(-W_disorder, W_disorder), nsite)
          # sz_arr = rand(Uniform(-W_disorder, W_disorder), nsite)
          # sx_arr = randn(nsite) .+ sinh(beta_h) * sqrt(n_alpha) * 1/4
          # sy_arr = randn(nsite)
          # sz_arr = randn(nsite)

          sx_arr = randn(nsite) * 1/2 .+ sinh(beta_h) * sqrt(n_alpha) * 1/4
          sy_arr = randn(nsite) * 1/2
          sz_arr = randn(nsite) * 1/2

          u0 = zeros(3*nsite)
          for j in 1:nsite
               # u0[3*(j-1)+1] = 1.0
               # u0[3*(j-1)+2] = 0.1*cos(j/nsite*2pi)
               # u0[3*(j-1)+3] = 0.1*sin(j/nsite*2pi)

               u0[3*(j-1)+1] = sx_arr[j]
               u0[3*(j-1)+2] = sy_arr[j]
               u0[3*(j-1)+3] = sz_arr[j]
          end

          
          # p[1:3] xi_1,2,3
          # p[4:nsite*(nsite+1)÷2+3] Jij
          # p[nsite*(nsite+1)÷2+4] nsite
          # p = zeros(nsite*(nsite+1)÷2+4)

          Jij_mat = zeros(nsite, nsite)

          index = 1
          for i in 1:nsite
               for j in i+1:nsite
                    temp_rand = (randn() * sigmasqrt)
                    Jij_mat[i,j] = temp_rand * sqrt(n_alpha)
                    index += 1
               end
          end
          # Jij_mat = 1/2 * (Jij_mat + transpose(Jij_mat))

          # No 1/2, since it's spin operator s_i = σ_i 1/2
          Jij_mat = (Jij_mat + transpose(Jij_mat))

          p = zeros(Jij_num+4)
          p[1:3] = xi_arr
          p[4:Jij_num+3] = vec(Jij_mat)
          p[Jij_num+4] = nsite

          prob = ODEProblem(parameterized_randomspin!,u0,tspan,p)
          sol = solve(prob,reltol=1e-5,adaptive=false,dt=delta_t)

          # plot(sol,vars=(1,2,3))
          # xs = (sol.t)*J

          t_len = length(sol.t)
          println("norm diff t_arr $(norm(sol.t - t_arr_fix))")
          u_data = sol.u
          Sxyz_t = zeros(t_len, 3)
          for i in 1:t_len
               for j in 1:nsite
                    Sxyz_t[i,1] += u_data[i][3*(j-1)+1] / nsite
                    Sxyz_t[i,2] += u_data[i][3*(j-1)+2] / nsite
                    Sxyz_t[i,3] += u_data[i][3*(j-1)+3] / nsite
               end
          end

          ##############################################################
          #2023/6/28
          # # Interpolations not works, since its better to use uniform data
          # times = (sol.t)*J
          # t = LinRange(0,1,length(times))
          # for a in 1:3
          #      vals = Sxyz_t[:,a]
          #      itp = Interpolations.scale(interpolate(hcat(times,vals), (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t, 1:2)
          #      ti = 0:1/400:1
          #      timesitp, valsitp = [itp(t,1) for t in ti], [itp(t,2) for t in ti]
          #      Sxyz_t_rand[:,a,irand] = valsitp

          #      display(timesitp)

          #      # interp = interp_cubic_Sxyz[a]
          #      # for i in 1:t_num_set
          #      #     Sxyz_t_rand[i,a,irand] = interp(t_arr_fix[i])
          #      # end
          # end
          #############################################################

          Sxyz_t_rand[:,:,irand] = Sxyz_t
     end


     label_Sxyz = Array{String, 2}(undef, 1, 3)
     label_Sxyz[1,1] = "\$ S^x \$"
     label_Sxyz[1,2] = "\$ S^y \$"
     label_Sxyz[1,3] = "\$ S^z \$"

     Sxyz_t_aver = mean(Sxyz_t_rand, dims=3)[:,:,1]

     fontsize = 11
     figplot = plot(t_arr_fix*J, Sxyz_t_aver, title = L"$ N=%$(nsite), J=%$(J), \xi=%$(xi_arr),\beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$tJ$",ylabel=L"\langle\hat{S}^{\alpha}(t)\rangle", label=label_Sxyz,  titlefontsize=fontsize)

     fileprefix = "xi$(xi_arr)_beta_h$(beta_h)_Nrand$(Nrealization)_tspan$(tspan)"

     # total_dirname = "figs_var1"
     # if !isdir(total_dirname)
     #      mkdir(total_dirname)
     # end

     savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix).pdf")
     savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix).png")
     JLD2.save("$(total_dirname)/Sxyz_meanfield_$(fileprefix).jld2", "tJ_arr", t_arr_fix*J, "Sxyz_t_aver",Sxyz_t_aver, 
          "Sxyz_t_rand",Sxyz_t_rand)

     return t_arr_fix*J, Sxyz_t_aver
end

# prob, sol, figplot = main()
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
# xi_arr_all = [[1.0, 0.8, 0.2]]


xi_arr_all = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 0.2], [-0.05, -0.15, 0.2], [0.0, -0.2, 0.2], [0.05, -0.25, 0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]]
# xi_arr_all = [[0.2, 0.2, -0.4], [0.15, 0.2, -0.35], [0.1, 0.2, -0.3], [0.05, 0.2, -0.25], [0.0, 0.2, -0.2], [-0.05, 0.2, -0.15], [-0.1, 0.2, -0.1], [-0.15, 0.2, -0.05], [-0.2, 0.2, 0.0], [-0.25, 0.2, 0.05], [-0.3, 0.2, 0.1], [-0.35, 0.2, 0.15], [-0.4, 0.2, 0.2]]

# xi_arr_all = [[0.2, -0.4, 0.2], [0.2, -0.35, 0.15], [0.2, -0.3, 0.1], [0.2, -0.25, 0.05], [0.2, -0.2, 0.0], [0.2, -0.15, -0.05], [0.2, -0.1, -0.1], [0.2, -0.05, -0.15], [0.2, 0.0, -0.2], [0.2, 0.05, -0.25], [0.2, 0.1, -0.3], [0.2, 0.15, -0.35], [0.2, 0.2, -0.4]]

J = 1.0
tspan = (0.0,16.0)
delta_t = 0.04
t_arr_fix = collect(tspan[1]:(delta_t):tspan[2])
t_num_set = length(t_arr_fix)
nsite = 2000
n_alpha = 16
# Jrescale = 5
# beta_h = 10*0.04 / Jrescale

beta_h = 10*0.04
Nrealization = 20

# total_dirname = "figs_exp_x_reduced_beta_h"
total_dirname = "figs_exp_x"
if !isdir(total_dirname)
     mkdir(total_dirname)
end

fontsize = 11

xi_len = length(xi_arr_all)
Sxyz_t_aver_xi_all = zeros(xi_len, t_num_set, 3)
for (i,xi_arr) in enumerate(xi_arr_all)
     tJ_arr_fix, Sxyz_t_aver = main(nsite, n_alpha, beta_h, Nrealization, J, xi_arr, tspan, delta_t, total_dirname)
     Sxyz_t_aver_xi_all[i,:,:] = Sxyz_t_aver
end

coef_fit_arr = zeros(length(xi_arr_all), 5)
fitstep_span = (1, 16*25)
for ll in 1:length(xi_arr_all)
    @. model(x, p) = p[1]*exp(-x*p[2])*cos(p[3]*x+p[4]) + p[5]
    xdata = (t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]
    ydata = (Sxyz_t_aver_xi_all[ll,:,1])[fitstep_span[1]:fitstep_span[2]]
    p0 = [0.543437,  1.95245,     3.58682, -0.954433, 0.0]
    lb = [0.2,0.0,0.0,-π,-0.1]
    ub = [2.0,4.0,6.0,π,0.1]
    fitfunc = curve_fit(model, xdata, ydata, p0, lower=lb, upper=ub)

    coef_fit_arr[ll,:] = coef(fitfunc)
end

fit_all_xi = zeros(length(xi_arr_all), (fitstep_span[2]-fitstep_span[1]+1))
for ll in 1:length(xi_arr_all)
    @. modelfit(x, p) = p[1]*exp(-x*p[2])*cos(p[3]*x+p[4])
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
plot(t_arr_fix*J, transpose(Sxyz_t_aver_xi_all[:,:,1]), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$tJ$",ylabel=L"\langle\hat{S}^{x}(t)\rangle", legend = false, titlefontsize=fontsize, legendfont=font(fontsize))
figplot = plot!(((t_arr_fix*J)[fitstep_span[1]:fitstep_span[2]]), transpose(fit_all_xi), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$tJ$",ylabel=L"\langle\hat{S}^{x}(t)\rangle", legend = false, titlefontsize=fontsize, legendfont=font(fontsize), linestyle=:dash)

fileprefix = "beta_h$(beta_h)_Nrand$(Nrealization)_tspan$(tspan)"

savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall.pdf")
savefig(figplot, "$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall.png")

println(coef_fit_arr[:,3])

Sxfreqtheory = [A_polynomial(xi) for xi in xi_arr_all]

# Sxfreqtheory = [3.52888, 1.73286, 3.04033, 2.69281, 2.08214, 2.9698, 1.55393, 3.16804, 1.85311, 1.87208, 1.22274, 0.634665, 2.07319, 1.14433, 1.19348, 0., 0., 0., 0., 0.]
# SxfreqKB = [3.34044, 1.89207, 2.72674, 2.75744, 1.9145, 2.66232, 1.4702, 3.09645, 
# 1.68609, 1.75667, 1.42522, 0.602516, 2.10143, 1.27365, 1.36014, 
# 0.0070625, 0.000278918, 0.000200615, 0.000639637, 0.350797]
scatter(hcat(coef_fit_arr[:,3], Sxfreqtheory), title = L"$ N=%$(nsite), J=%$(J), \beta h=%$(beta_h),N_{\mathrm{rand}}=%$(Nrealization)$", xlabel=L"$\mathrm{Cases}$",ylabel=L"\Omega_x", label=["meanfield" "\$J\\sqrt{A}\$"],  titlefontsize=fontsize, legendfont=font(fontsize))
savefig("$(total_dirname)/Sxfreq_meanfield_$(fileprefix)_xiall.pdf")
savefig("$(total_dirname)/Sxfreq_meanfield_$(fileprefix)_xiall.png")


JLD2.save("$(total_dirname)/Sxyz_meanfield_$(fileprefix)_xiall.jld2", "tJ_arr", t_arr_fix*J, "Sxyz_t_aver_xi_all", Sxyz_t_aver_xi_all, "coef_fit_arr", coef_fit_arr)