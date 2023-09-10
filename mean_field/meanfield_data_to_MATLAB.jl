using Plots
using JLD2
using MAT

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


xi_arr_all = xi_arr_all_x
adirect = 1

# scaling_c = [0.951, 0.948, 0.899]
scaling_c = [0.9876, 0.9869, 0.9054]

xi_len = length(xi_arr_all)
Sxyz_t_aver_xi_all = zeros(xi_len, t_num_set, 3)
total_dirname = "figs_exp_$(direct_string[adirect])"
# # fitstep_span = (1, 12*25)
# # fitstep_span = (2, 12*25)

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

Sx_t_aver_xi_all = zeros(xi_len, t_num_set)
Sx_t_aver_xi_all = Sxyz_t_aver_xi_all[:,:,1]
for l in 1:13
    Sx_t_aver_xi_all[l,:] = Sx_t_aver_xi_all[l,:] / Sx_t_aver_xi_all[l,1]
end


file = matopen("MeanfieldSx.mat", "w")
write(file, "Expression1", Sx_t_aver_xi_all)
close(file)


file = matopen("MeanfieldTime.mat", "w")
write(file, "Expression1", t_arr_fix)
close(file)