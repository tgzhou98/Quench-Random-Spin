using Distributed
using JLD2
using DelimitedFiles
using Plots
using LaTeXStrings
# gr()
# using Arpack
ENV["GKSwstype"]="100"
gr()

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


############################# Dynamics


function disorder_aver_Sx(disorder_num::Int64, beta::Float64, N::Int64, J::Float64, Jbar::Float64,
     f::Array{Float64,1}, hinJ::Float64, quench_hinJ::Float64, t_list::Array{Float64, 1})
    # disorder_num = 100
    Sx_list_all_disorder  = zeros(length(t_list), disorder_num)

    println("f=$(f), beta=$(beta), N=$(N), Jbar=$(Jbar), hi=$(hinJ), hf=$(quench_hinJ)")
    # for i in 1:disorder_num
    # @showprogress pmap(1:disorder_num) do i
    for i in 1:disorder_num
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

function main()

    N = 5
    # disorder_num = 10
    disorder_num = 1000
    # 
    # N = 6
    # disorder_num = 50
    # 
    # t_list only use the Delta_t


    data_dirname = "rs_Sx_collapse_dynamics_smallSx_20220707"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end

    # hinJ = 10.0
    # quench_hinJ = 0.0
    # h = 10.0
    h = 0.1
    quench_h = 0.0


    xi_len = 2*2 + 1
    xi_step = 0.5

    Jlist = zeros(xi_len)
    flist = zeros(xi_len,3)
    # ij_map = zeros(Int64, xi_len*xi_len, 2)
    # initial flist_len
    flist_len = 1
    for i in 1:xi_len
        flist[flist_len, 1] = 1
        flist[flist_len, 2] = (i - (xi_len - 1)÷2 -1 + 1/2) * xi_step
        flist[flist_len, 3] = - 1/(2*flist[flist_len, 2])
        
        # ij_map[flist_len, 1] = i
        # ij_map[flist_len, 2] = j
        Jlist[flist_len] = 1 / norm(flist[flist_len,:])

        flist_len += 1
    end
    flist_len -= 1

    # Deltalist = [0.0, -2.0, 2.0]
    # flist = [[1.0, -0.5, 0.0], [1.0, -0.8, -0.5], [1.0, -1.0, -1.0], [1.0, 1.0, 0.5], [1.0, 1.0, -0.5], [1, 1, -2], [1, 1, -1], [1, 1, 5], [1, 1, 10], [1, 1, -1.5], [1, 1, -2.5], [1, 1, -3.0], [1.0, -2.0, 0.0], [1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]
    # Deltalist = [-0.5, 0.0, 0.5]
    # betalist = [0.0, 4.0, 15.0, 20.0]
    # betalist = [4.0, 30.0]
    betalist = [0.1]
    # shorttime
    # t_list = collect(0:0.05:50)
    # longtime
    # t_list = collect(0:0.5:30)


    tJ_list = collect(0:0.02:20)

    Sx_aver_listall = zeros(length(tJ_list), flist_len, length(betalist))

    t_listall =  zeros(length(tJ_list), flist_len)

    Jbar = 0.0
    for i in 1:flist_len
    # for Jbar in [-0.3]
        f = flist[i,:]
        J = Jlist[i]
        t_list = tJ_list / J
        t_listall[:,i] = t_list

        hinJ = h / J
        quench_hinJ = quench_h / J
        for (j, beta) in enumerate(betalist)
            # Delta = -0.73
            # Jbar = 0.0
            Sx_aver_list, Sx_list_all_disorder = disorder_aver_Sx(disorder_num, beta, N, J, Jbar, f, hinJ, quench_hinJ, t_list)
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
        JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall, "t_list", t_list)
    end

    writedlm("$(data_dirname)/f_arr.csv", flist[:,:], ",")
    writedlm("$(data_dirname)/Sx_aver.csv", Sx_aver_listall[:,:,1], ",")
    writedlm("$(data_dirname)/tJ_list.csv", tJ_list, ",")

    # JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall, "t_list", tperbeta_list*betalist[1])
    # 
    # 
    # Sx_aver_listall = JLD2.load("Sx_aver.jld2", "Sx_aver")


    #################################################################################
    # beta as legends, f as two figures
    label_f = Array{String, 2}(undef, 1, flist_len)
    for ll in 1:flist_len
        label_f[1,ll] = "\$ f=$(flist[ll,:]) \$"
    end


    maxplotstep = length(tJ_list)
    # p1 = plot(tJ_list[1:maxplotstep], Sx_aver_listall[1:maxplotstep,:,1], title = "\$ \\beta=$(betalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(h), h_f=$(quench_h) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
    p1 = plot(t_listall[1:maxplotstep,:], Sx_aver_listall[1:maxplotstep,:,1], title = "\$ \\beta=$(betalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(h), h_f=$(quench_h) \$", xlabel="\$ t \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
    # p2 = plot(tperbeta_list[1:maxplotstep]*betalist[1], Sx_aver_listall[1:maxplotstep,:,2], title = "\$ \\beta=$(betalist[2]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
    # plot(p1, p2, layout=(2))
    plot(p1)

    # JLD2.save("$(data_dirname)/Sx_aver.jld2", "Sx_aver", Sx_aver_listall)
    savefig("$(data_dirname)/thermalSx.pdf")

    ###################################### END Dynamics ###############################

    
    return tJ_list, Sx_aver_listall, flist
end

tJ_list, Sx_aver_listall, flist = main()