include("run_Jbar.jl")

######################################################################################################
# ONLY for debug
# Single shot run in global scope
######################################################################################################
########################## Setting 1
# ρω, abscissaeω, Sxyz_new, h0, Jbar=equilibrium_Jbar(4.0, 0.05, "x", -0.73, 1.0, 2π, 8.0, 50.0)
# J, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
#     dynamics_init_Jbar(2.0, "x", -0.73, 2π, ρω, 8.0, 50.0, 1001);
########################## Setting 2
# init_GRω_path = "./ρω.jld"
# GRω = JLD.load(init_GRω_path, "GRω")
# ρω, abscissaeω, Sxyz_new, h0, Jbar=equilibrium_Jbar(15.0, 0.2, "x", 0.0, 1.0, 15.0/5.0*2π, 2.0, 15*3.0, GRω);
# J, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
#     dynamics_init_Jbar(15.0, "x", 0.0, 15.0/5.0*2π, ρω, 2.0, 15*3.0, 1001);
# baym_kadanoff_Jbar(J, h0, "x",Sxyz_new, Jbar, nt, Δt, f
#     ,GGreaterDynamicsttp
#     ,GLesserDynamicsttp
#     ,ΣGreaterDynamicsttp
#     ,ΣLesserDynamicsttp
#     ,GRetardedDynamicsttp
#     ,GAdavancedDynamicsttp
#     ,ΣRetardedDynamicsttp
#     ,ΣAdavancedDynamicsttp
#     ,"TI"
#     ,false
# )
# plot_Sx(GLesserDynamicsttp,   2.0, 0.05, "x", -0.73, 800, 1250, 50/2000, 1001)
# plot_Gt(GLesserDynamicsttp, 2.0, 0.05, "x", -0.73, 800, 1250, 50/2000, 1001)

# Debug dynamics

# get_dynmaics_Jbar(4.0, 0.2, "x", -0.73, 0.0, 2π,
#     8.0, 50.0, 1001, "quench")

# init_GRω_path = "initial_guess_fixβ_J_tune_init_GRω_finecutoff_MQ_march30_wormhole/Jbar-1_0/betaJ27_0_hinJ0_2/ρω.jld"
# GRω = JLD.load(init_GRω_path, "GRω")

# get_dynmaics_Jbar(29.0, 0.2, "x", 0.0, 0.0, 29.0, 3.0, 7.0*29, 4001, "TI", GRω)
# mv("quench_data_Jbar", "get_dynmaics_Jbar_alpha_MQ_Jbar0_TI_april25_test")

############################ April 4, test imaginary spectrum
# Gt_init = JLD.load("./imagtime_Jbar_init_guess/Jbar_L/imagtime_eq.jld", "Gt")
# # abscissaet, Sxyz_new, h0, free_energy, Jbar, Gt=equilibrium_Jbar_imagtime(15.0, 0.2, "x", -0.73, 1.0, 15.0, Gt_init);
# abscissaet, Sxyz_new, h0, free_energy, Jbar, Gt=equilibrium_Jbar_imagtime(4.0, 0.0, "x", -0.73, 1.0, 4.0);

############################ may15, test Sx iteration, h=0

# init_GRω_path = "./ρω.jld"
# GRω = JLD.load(init_GRω_path, "GRω")
# # GRω = Array{ComplexF64, 3}(undef,0,0,0)
# ρω, abscissaeω, Sxyz_new, h0, Jbar=equilibrium_Jbar(12.0, 0.2, "x", -0.73, -0.3, 2π, 8.0, 50.0, GRω, true)