include("run_Jbar.jl")

####################### dynamics_20220428_hightemp_SpinModel_Jbar0 ###################
# For various T, h and f 
# Search all phase diagram
########################################


# init_GRω_path = "spectrum_anal_Jbar_init_guess/betaJ13_0_hinJ0_1/ρω.jld"
# init_GRω_path = "ρω.jld"
init_GRω_path = "dynamics_20230715_draftRandomf_hightemp_SpinModel_Jbar0_rhox/set_5/Jbar0_0/ToverJ25_0_hinJ10_0/ρω.jld"
# init_GRω_path = "phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff_MQ_march30_wormhole_fortest/Jbar-0_5/betaJ31_0_hinJ0_2/ρω.jld"

#init_GRω_path = "initial_guess_fixβ_J_tune_init_GRω_finecutoff_MQ_march30_wormhole/Jbar-1_0/betaJ15_0_hinJ0_2/ρω.jld"

GRω = JLD.load(init_GRω_path, "GRω")


# total_dirname = "dynamics_Jbar_CutoffNum4001_20220412_fixbeta_h_Jbar0_finetune"
# total_dirname = "dynamics_20220622_hightemp_SpinModel_Jbar0_anyf"
# total_dirname = "dynamics_20221118_hightemp_SpinModel_Jbar0_anyf"
# total_dirname = "dynamics_20230313_draftRandomf_hightemp_SpinModel_Jbar0_anyf"
total_dirname = "dynamics_20230715_draftRandomf_hightemp_SpinModel_Jbar0_rhox"
if !isdir(total_dirname)
    mkdir(total_dirname)
end


###############################################
# low T
# hlist = [0.016, 0.03, 0.05, 0.05]
# Λt_arrlist = [100.0, 100.0, 100.0, 100.0]
# Λω = 4.0


# hlist = [0.05]
# Λt_arrlist = [100.0]
# Λω = 5.0
###############################################

# 2023/03/13
#=
flist = [[-1.46, 1.89, -1.25], [1.47, -1.78, 0.27], [0.36, -1.66, 
  0.88], [1.79, -1.68, 0.99], [-0.8, 1.07, -0.76], [-0.51, 
  0.97, -1.51], [0.25, 0.08, -1.42], [-1.58, 
  1.05, -1.88], [-0.01, -0.29, 1.34], [0.09, -1.54, 
  0.18], [1.51, -0.22, 1.54], [0.25, 0.06, 0.81], [-1.49, -1.27, 
  0.83], [0.89, -0.1, -1.66], [-0.73, 1.89, 0.22], [0.56, 0.58, 
  0.93], [-1.83, -0.97, -1.04], [-1.15, 1.23, 
  1.27], [1.31, -1.68, -1.31], [1.2, 1.07, 0.]]
=#
# flist = [[1.0, 0.8, -1.5], [1.0, 0.8, 1.5]]
# flist = [[1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]
# flist = [[1.0, 1.0, -2.0], [1.0, -1.0, 0.0], [1.0, 1.0, 2.0]]
# flist = [[1.0, 1.0, -1.0], [1.0, 0.5, -0.8], [1.0, -0.5, 0.8]]
# flist = [[1.0, -0.5, 0.0], [1.0, -1.0, 0.0], [1.0, -2.0, 0.0], [1.0, -1.0, -1.0]] # dynamics_20220619_lowtemp_SpinModel_Jbar0


# 2023/07/15
Jrescale = 5
flist = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 
   0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 
   0.2], [-0.05, -0.15, 0.2], [0.02, -0.22, 0.2], [0.05, -0.25, 
   0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]] * Jrescale
# flist = [[0.02, -0.22, 0.2]] * 5 # is not converge
# flist = [[0.05, -0.25, 
#    0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]] * Jrescale

###############################################
# high T
# hlist = [10.0, 10.0, 10.0, 10.0]
# Λt_arrlist = [20.0, 20.0, 20.0, 20.0]
# Λω = 12.0
hlist = [10.0 for i in 1:length(flist)]
Λt_arrlist = [20.0 for i in 1:length(flist)]
Λω = 12.0  # Standard cutoff
# Λω = 8.0
###############################################

# hlist = [0.0, 10.0]
# Tlist = [10.0/0.4, 10.0/0.4]

Tlist = hlist / 0.4


# for f_ind in 1:2
# flist = [[1.0, 1.0, -2.0]]
# for set_index in 1:length(hlist)

writedlm("$(total_dirname)/flist.csv", flist, ',')
writedlm("$(total_dirname)/T_h_list.csv", hcat(Tlist, hlist), ',')
writedlm("$(total_dirname)/cutofft_arrlist.csv", Λt_arrlist, ',')

for set_index in 1:length(hlist)
    set_dirname = "set_$(set_index)"
    # data_dirname = "H_1"
    if !isdir(set_dirname)
        mkdir(set_dirname)
    end

    f = flist[set_index]




    #ξ

    # begin_betaJ = 8.0
    # end_betaJ = 20.0
    # step_betaJ = 4.0

    # begin_betaJ = 5.0
    # end_betaJ = 5.0
    # step_betaJ = 6.25
    # begin_hinJ = 1.0
    # end_hinJ = 1.0
    begin_hinJ = hlist[set_index]
    end_hinJ = hlist[set_index]
    step_hinJ = 0.4
    field = "x"
    # J is 0.75
    std_betaJ = 2π 

    # alpha = -0.5
    cutoffNum = 4001
    # cutoffNum = 501
    # cutoffNum = 4501
    # Tlist_alpha_list = [0.2]
    # Tlist_alpha_list = [0.2, 0.5, 2.0]
    # Tlist_alpha_list = [1/(0.005)]
    # Tlist_alpha_list = [1/(0.00625)]
    Tlist_alpha_list = [Tlist[set_index]]
    betaJ_num = length(Tlist_alpha_list)
    betaJlist = 1.0 ./ Tlist_alpha_list
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ
    Λt_arr = ones(betaJ_num, h_num)

    for h_ind in 1:h_num
        # Λt_arr[:,h_ind] = [2500.0 * betaJ  for betaJ in betaJlist]
        Λt_arr[:,h_ind] = [Λt_arrlist[set_index]]
    end

    #Λt_arr[betaJ_num-1,1] = 6.0 * betaJlist[betaJ_num-1]
    #Λt_arr[betaJ_num,1] = 7.0 * betaJlist[betaJ_num]


    #Λt_arr[betaJ_num-3,1] = 6.0 * betaJlist[betaJ_num-3]
    #Λt_arr[betaJ_num-2,1] = 7.0 * betaJlist[betaJ_num-2]
    #Λt_arr[betaJ_num-1,1] = 7.0 * betaJlist[betaJ_num-1]
    #Λt_arr[betaJ_num,1] = 9.0 * betaJlist[betaJ_num]


    #for Jbar_val in [1.0, 0.0]
    for Jbar_val in [0.0]
        Jbar_str = replace("$(Jbar_val)", "."=>"_")
        run_phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff(
            Tlist_alpha_list,begin_hinJ ,end_hinJ ,step_hinJ, std_betaJ, field, f, Jbar_val, Λω, Λt_arr, cutoffNum, GRω, "quench", false)

        # phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff_generate(
        #        begin_betaJ ,end_betaJ ,step_betaJ ,begin_hinJ ,end_hinJ ,step_hinJ, std_betaJ, field, alpha, Jbar_val, Λω, Λt_arr, 4001, GRω, "quench", false)

        mv("phase_diagram_dynamics", "$(set_dirname)/Jbar$(Jbar_str)", force=true)
    end


    mv("$(set_dirname)", "$(total_dirname)/$(set_dirname)", force=true)
end
