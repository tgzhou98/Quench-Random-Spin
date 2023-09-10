using Distributed
using JLD2
using DelimitedFiles
using Plots
# gr()
# using Arpack
ENV["GKSwstype"]="100"
gr()

@everywhere begin
    # using ProgressMeter
    using Distributed
    using SharedArrays

    include("rd_model_ed.jl")

    function one_shot(beta::Float64, N::Int64, J::Float64, Jbar::Float64, f::Array{Float64,1}, hinJ::Float64, quench_hinJ::Float64, t_list::Array{Float64, 1})
        # beta = 15
        # N = 8
        # println(f)
        # Jij_term_ops_arr = generate_Jij_term_ops_arr(N, f)
        # HamJij = RD_ham_Jij(N, J, Jbar, f, Jij_term_ops_arr)
        HamJij = generate_RD_ham_Jij(N, J, Jbar, f)
        Hamh_bare = RD_ham_h(N, 1.0, 1.0)
        Hamh0 = Hamh_bare * J * hinJ
        Ham1 = HamJij + Hamh0
        rho = exp(- beta*Ham1)
        rho = rho / tr(rho)
        # sqrtrho = sqrt(rho)

        a = 1
        Sx_op = zeros(ComplexF64, 2^N, 2^N)
        # Sx_op_list = zeros(ComplexF64, 2^N, 2^N, N)
        for i in 1:N
            hi_term_ops = fill(Idmatrix, N)
            hi_term_ops[i] = PauliMatrices[a]
            ####################################
            # Bug, not use 1/2 here
            ####################################
            # Sx_op_list[:,:,i] = 1 / 2 * foldl(⊗, hi_term_ops)
            Sx_op[:,:] += foldl(⊗, hi_term_ops) / 2
        end

        # initial equilibrium hamiltonian
        Ham_quench = HamJij + Hamh_bare * quench_hinJ * J
        # rho_t = rho

        Sx_list = zeros(length(t_list))
        deltaT = t_list[2] - t_list[1]

        # prepare a bare h hamiltonian

        # DON'T USE MOLECULAR FIELD APPRXIMATION
        vals, vecs = eigen(Ham_quench)
        for (i, t) in enumerate(t_list)
            # vals, vecs = eigen(Ham_quench)

            expEigenvalsH1 = [exp(-im * t * ev) for ev in vals]
            DH1 = Diagonal(expEigenvalsH1)

            U_H1t = vecs * DH1 #exp(-itD) *
            U_H1t = BLAS.gemm('N','C', U_H1t, vecs)

            # ############### DON'T USE THIS MOLECULAR FILED CODE
            # U_H1t = exp(-im * deltaT * Ham_quench)
            # ###################################################

            U_H2t = adjoint(U_H1t) #exp(itD) *
            Sxt = U_H2t*Sx_op*U_H1t


            Sx_list[i] = real(tr(rho*Sxt)) / N

            # if i == 1
            #     println(Sx_list[i])
            # end

            # ############### DON'T USE THIS MOLECULAR FILED CODE
            # # Update Ham_quench
            # Hamh_new = RD_ham_h(N, J, -Sx_list[i]*Jbar*J)
            # Hamh_new = Hamh_bare * (-Sx_list[i]*Jbar*J) 
            # Ham_quench = HamJij + Hamh_new
            # ###################################################
        end

        return Sx_list
    end

end

############################# Dynamics


function disorder_aver_Sx(disorder_num::Int64, beta::Float64, N::Int64, J::Float64, Jbar::Float64,
     f::Array{Float64,1}, hinJ::Float64, quench_hinJ::Float64, t_list::Array{Float64, 1})
    # disorder_num = 100
    Sx_list_all_disorder  = SharedArray{Float64}(length(t_list), disorder_num)

    println("f=$(f), beta=$(beta), N=$(N), Jbar=$(Jbar), hi=$(hinJ), hf=$(quench_hinJ)")
    # for i in 1:disorder_num
    @showprogress pmap(1:disorder_num) do i
        if i % 10 == 1
            println("Disorder step: $(i)")
        end
        Sx_list_all_disorder[:, i] = one_shot(beta, N, J, Jbar, f, hinJ, quench_hinJ, t_list)
    end
    println("finished\n\n")
    
    Sx_list_final = mean(Sx_list_all_disorder, dims=2)
    return Sx_list_final, Sx_list_all_disorder
end


# N = 10
# J = 1
# Delta = -0.73
# hinJ = 0.2
# HamJij = RD_ham_Jij(N, J, Delta)
# Hamh0 = RD_ham_h(N, J, hinJ)
# Ham1 = HamJij + Hamh0

# # Slower than default exp function
# function tempfunc(Ham_quench, deltaT)
#     vals, vecs = eigen(Ham_quench)

#     expEigenvalsH1 = [exp(-im * deltaT * ev) for ev in vals]
#     DH1 = Diagonal(expEigenvalsH1)

#     U_H1t = vecs * DH1 #exp(-itD) *
#     # U_H1t = BLAS.gemm('N','C', U_H1t, vecs)
#     U_H1t = exp(-im * deltaT * Ham_quench)
# end
# @btime exp(-1im*1*Ham1)

###################################### BEGIN Dynamics ###############################
# t_list = collect(0:0.5:30)
# # for Jbar in [0.0]
# for Jbar in [0.0, 1.0, -0.3]
#     for beta in [15.0, 31.0]
#     # for beta in [31.0, 36.0]
#         Delta = 0.0
#         # Jbar = 0.0
#         Sx_aver_list, Sx_list_all_disorder = disorder_aver_Sx(disorder_num, beta, N, 1.0, Jbar, Delta, hinJ, quench_hinJ, t_list)
#         # plot(t_list, Sx_list_all_disorder) Δ
#         fileprefix = "SxDelta$(Delta)_beta$(beta)_N$(N)_Jbar$(Jbar)"

#         plot(t_list, Sx_aver_list, title = "\\Delta =$(Delta), beta=$(beta), N=$(N), Jbar=$(Jbar)", xlabel="\$ tJ \$",ylabel="\$ S_x \$")
#         savefig("$(data_dirname)/$(fileprefix).pdf")

#         plot(t_list, Sx_list_all_disorder, title = "\\Delta =$(Delta), beta=$(beta), N=$(N), Jbar=$(Jbar)", xlabel="\$ tJ \$",ylabel="\$ S_x \$")
#         savefig("$(data_dirname)/$(fileprefix)_alldisorder.pdf")
#         JLD2.save("$(data_dirname)/$(fileprefix).jld2", "Sx_list_all_disorder", Sx_list_all_disorder.s)
#         writedlm("$(data_dirname)/$(fileprefix).csv", hcat(t_list, Sx_aver_list) , ',')
#     end
# end


N = 8
# disorder_num = 10
disorder_num = 1000
# 
# N = 6
# disorder_num = 50
# 
# t_list only use the Delta_t


# data_dirname = "rs_Sx_dynamics_20221121"
# data_dirname = "rs_Sx_draftRandomf_dynamics_20230314"
data_dirname = "rs_Sx_draftRandomf_dynamics_20230716_rhox"
if !isdir(data_dirname)
    mkdir(data_dirname)
end

hinJ = 10.0
quench_hinJ = 0.0


# Deltalist = [0.0, -2.0, 2.0]
# flist = [[1.0, -0.5, 0.0], [1.0, -0.8, -0.5], [1.0, -1.0, -1.0], [1.0, 1.0, 0.5], [1.0, 1.0, -0.5], [1, 1, -2], [1, 1, -1], [1, 1, 5], [1, 1, 10], [1, 1, -1.5], [1, 1, -2.5], [1, 1, -3.0], [1.0, -2.0, 0.0], [1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]
# flist = [[1.0, 0.8, -1.5], [1.0, 0.8, 1.5], [1.0, 1.0, -2.0], [1.0, 1.0, 2.0]]
# flist = [[1.0, 1.0, -0.5], [1.0, 1.0, -1.5], [1.0, 1.0, -2.5], [1.0, 1.0, -3.0]]
# flist = [[-1.46, 1.89, -1.25], [1.47, -1.78, 0.27], [0.36, -1.66, 
#   0.88], [1.79, -1.68, 0.99], [-0.8, 1.07, -0.76], [-0.51, 
#   0.97, -1.51], [0.25, 0.08, -1.42], [-1.58, 
#   1.05, -1.88], [-0.01, -0.29, 1.34], [0.09, -1.54, 
#   0.18], [1.51, -0.22, 1.54], [0.25, 0.06, 0.81], [-1.49, -1.27, 
#   0.83], [0.89, -0.1, -1.66], [-0.73, 1.89, 0.22], [0.56, 0.58, 
#   0.93], [-1.83, -0.97, -1.04], [-1.15, 1.23, 
#   1.27], [1.31, -1.68, -1.31], [1.2, 1.07, 0.]]


# 2023/7/16
Jrescale = 5
flist = [[-0.4, 0.2, 0.2], [-0.35, 0.15, 0.2], [-0.3, 0.1, 0.2], [-0.25, 
   0.05, 0.2], [-0.2, 0.0, 0.2], [-0.15, -0.05, 0.2], [-0.1, -0.1, 
   0.2], [-0.05, -0.15, 0.2], [0.0, -0.2, 0.2], [0.05, -0.25, 
   0.2], [0.1, -0.3, 0.2], [0.15, -0.35, 0.2], [0.2, -0.4, 0.2]] * Jrescale

# flist = [[1.0, 0.8, -1.5], [1.0, 0.8, 1.5]]
# Deltalist = [-0.5, 0.0, 0.5]
# betalist = [0.0, 4.0, 15.0, 20.0]
# betalist = [4.0, 30.0]
betalist = [1/(hinJ/0.4)]
# shorttime
# t_list = collect(0:0.05:50)
# longtime
# t_list = collect(0:0.5:30)


tperbeta_list = collect(0:0.01:12) / betalist[1]

Sx_aver_listall = zeros(length(tperbeta_list), length(flist), length(betalist))

writedlm("$(data_dirname)/f_arr.csv", flist[:,:], ",")
writedlm("$(data_dirname)/t_list.csv", tperbeta_list*betalist[1], ",")

Jbar = 0.0
for (i, f) in enumerate(flist)
# for Jbar in [-0.3]
    for (j, beta) in enumerate(betalist)
        t_list = tperbeta_list * beta
        # Delta = -0.73
        # Jbar = 0.0
        Sx_aver_list, Sx_list_all_disorder = disorder_aver_Sx(disorder_num, beta, N, 1.0, Jbar, f, hinJ, quench_hinJ, t_list)
        Sx_aver_listall[:,i,j] = Sx_aver_list
        # plot(t_list, Sx_list_all_disorder) α
        # fileprefix = "SxDelta$(Delta)_beta$(beta)_N$(N)_Jbar$(Jbar)_hi$(hinJ)_hf$(quench_hinJ)"

        # println(Sx_aver_list)

        # plot(t_list, Sx_aver_list, title = "\$ \\Delta =$(Delta), \\beta=$(beta), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$" ,xlabel="\$ tJ \$",ylabel="\$ \\mathrm{G} \$")
        # savefig("$(data_dirname)/$(fileprefix).pdf")

        # plot(t_list, Sx_list_all_disorder, title = "\$ \\Delta =$(Delta), \\beta=$(beta), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel="\$ \\mathrm{G} \$")
        # savefig("$(data_dirname)/$(fileprefix)_alldisorder.pdf")
        # JLD2.save("$(data_dirname)/$(fileprefix).jld2", "Sx_list_all_disorder", Sx_list_all_disorder.s)
        # writedlm("$(data_dirname)/$(fileprefix).csv", hcat(t_list, Sx_aver_list) , ',')
    end
    JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall, "t_list", tperbeta_list*betalist[1])
    writedlm("$(data_dirname)/Sx_aver.csv", Sx_aver_listall[:,:,1], ",")
end


# JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall, "t_list", tperbeta_list*betalist[1])
# 
# 
# Sx_aver_listall = JLD2.load("Sx_aver.jld2", "Sx_aver")


#=
#################################################################################
# beta as legends, Delta as two figures
label_beta = Array{String, 2}(undef, 1, length(betalist))
for ll in 1:length(betalist)
    label_beta[1,ll] = "\$ \\beta=$(betalist[ll]) \$"
end
using LaTeXStrings
# p1 = plot(tperbeta_list, Sx_aver_listall[:,1,:], title = "\$ \\Delta=$(Deltalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# p2 = plot(tperbeta_list, Sx_aver_listall[:,2,:], title = "\$ \\Delta=$(Deltalist[2]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# p3 = plot(tperbeta_list, Sx_aver_listall[:,3,:], title = "\$ \\Delta=$(Deltalist[3]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# p4 = plot(tperbeta_list, Sx_aver_listall[:,4,:], title = "\$ \\Delta=$(Deltalist[4]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
maxplotstep=100
p1 = plot(tperbeta_list[1:maxplotstep], Sx_aver_listall[1:maxplotstep,1,:], title = "\$ \\Delta=$(Deltalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
p2 = plot(tperbeta_list, Sx_aver_listall[:,2,:], title = "\$ \\Delta=$(Deltalist[2]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
p3 = plot(tperbeta_list, Sx_aver_listall[:,3,:], title = "\$ \\Delta=$(Deltalist[3]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
p4 = plot(tperbeta_list[1:maxplotstep], Sx_aver_listall[1:maxplotstep,4,:], title = "\$ \\Delta=$(Deltalist[4]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ t/\\beta \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# p4 = plot(legend=false,grid=false,foreground_color_subplot=:white)
plot(p1,p2,p3,p4, layout=(2,2))

# JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall)
savefig("$(data_dirname)/thermalSx.pdf")
# savefig("thermalSx.pdf")

# Sx_aver_listall = JLD2.load("$(data_dirname)/Sx_aver.jld2", "Sx_aver")
# Sx_aver_listall_normalize = 0.0*Sx_aver_listall
# for (i, t_list) in enumerate(t_list)
#     for (j, beta) in enumerate(betalist)
#     # for Jbar in [-0.3]
#         for (k, Delta) in enumerate(Deltalist)
#             Sx_aver_listall_normalize[i,j,k] = Sx_aver_listall[i,j,k] / Sx_aver_listall[1,j,k]
#         end
#     end
# end
# p1 = plot(t_list, Sx_aver_listall_normalize[:,1,:], title = "\$ \\Delta=$(Deltalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
# p2 = plot(t_list, Sx_aver_listall_normalize[:,2,:], title = "\$ \\Delta=$(Deltalist[2]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
# p3 = plot(t_list, Sx_aver_listall_normalize[:,3,:], title = "\$ \\Delta=$(Deltalist[3]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x_1(t) \hat{S}^x_1(0)) ", label=label_beta, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
# p4 = plot(legend=false,grid=false,foreground_color_subplot=:white)
# plot(p1,p2,p3,p4, layout=(2,2))
# savefig("$(data_dirname)/thermalSx_normalize.pdf")
#################################################################################
=#



#################################################################################
# beta as legends, f as two figures
label_f = Array{String, 2}(undef, 1, length(flist))
for ll in 1:length(flist)
    label_f[1,ll] = "\$ f=$(flist[ll]) \$"
end
using LaTeXStrings


maxplotstep=length(tperbeta_list)
p1 = plot(tperbeta_list[1:maxplotstep]*betalist[1], Sx_aver_listall[1:maxplotstep,:,1], title = "\$ \\beta=$(betalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# p2 = plot(tperbeta_list[1:maxplotstep]*betalist[1], Sx_aver_listall[1:maxplotstep,:,2], title = "\$ \\beta=$(betalist[2]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
# plot(p1, p2, layout=(2))
plot(p1)

# JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall)
savefig("$(data_dirname)/thermalSx.pdf")

###################################### END Dynamics ###############################
