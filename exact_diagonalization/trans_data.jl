using JLD
using DelimitedFiles


data_dirname = "rs_Sx_dynamics_20220622"

Sx_aver = JLD.load("$(data_dirname)/Sx_aver.jld", "Sx_aver")
t_list = JLD.load("$(data_dirname)/Sx_aver.jld", "t_list")

flist = [[1.0, -0.5, 0.0], [1.0, -0.8, -0.5], [1.0, -1.0, -1.0], [1.0, 1.0, 0.5], [1.0, 1.0, -0.5], [1, 1, -2], [1, 1, -1], [1, 1, 5], [1, 1, 10], [1, 1, -1.5], [1, 1, -2.5], [1, 1, -3.0], [1.0, -2.0, 0.0], [1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]

writedlm("$(data_dirname)/Sx_aver.csv", Sx_aver[:,:], ",")
writedlm("$(data_dirname)/t_list.csv", t_list, ",")
writedlm("$(data_dirname)/f_arr.csv", flist, ",")
