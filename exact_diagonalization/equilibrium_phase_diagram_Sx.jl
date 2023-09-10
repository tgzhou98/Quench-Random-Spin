using Distributed
using JLD
using DelimitedFiles
using Plots
# gr()
# using Arpack
ENV["GKSwstype"]="100"
gr()

include("rd_model_ed.jl")

function one_shot_Sx(beta::Float64, N::Int64, J::Float64, Jbar::Float64, f::Array{Float64,1}, hinJ::Float64)
    # beta = 15
    # N = 8
    # println(f)
    # Jij_term_ops_arr = generate_Jij_term_ops_arr(N, f)
    # HamJij = RD_ham_Jij(N, J, Jbar, f, Jij_term_ops_arr)
    HamJij = generate_RD_ham_Jij(N, J, Jbar, f)
    Hamh_bare = RD_ham_h(N, 1.0, 1.0)
    Hamh0 = Hamh_bare * J * hinJ
    Ham1 = HamJij + Hamh0
    # sqrtrho = sqrt(rho)
    rho = exp(- beta*Ham1)
    rho = rho / tr(rho)

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

    Sx_eq = real(tr(rho*Sx_op)) 

    return Sx_eq
end

############################# Dynamics


function random_realization_Sx(beta::Float64, N::Int64, J::Float64, Jbar::Float64,  hinJ::Float64, nrand_num::Int64, flist_arr::Array{Float64,2}, f_num::Int64)
    sigmasqrt = sqrt(4*J^2/N)
    H_XXYYZZ = generate_RD_ham_Jij_pauliXXYYZZ(N, J, Jbar)

    Hamh_bare = RD_ham_h(N, 1.0, 1.0)
    Hamh0 = Hamh_bare * J * hinJ

    # Sx operator
    Sx_op = zeros(ComplexF64, 2^N, 2^N)
    # Sx_op_list = zeros(ComplexF64, 2^N, 2^N, N)
    a = 1
    for i in 1:N
        hi_term_ops = fill(Idmatrix, N)
        hi_term_ops[i] = PauliMatrices[a]
        ####################################
        # Bug, not use 1/2 here
        ####################################
        # Sx_op_list[:,:,i] = 1 / 2 * foldl(⊗, hi_term_ops)
        Sx_op[:,:] += foldl(⊗, hi_term_ops) / 2
    end

    
    Sx_flist_val = zeros(f_num)
    for i in 1:f_num
        println("N=$(N), J=$J, f=$(flist_arr[i,:]),  step=$i of $(f_num)")
        HamJ_noRandJ = flist_arr[i,1]*H_XXYYZZ[:,:,1] + flist_arr[i,2]*H_XXYYZZ[:,:,2] + flist_arr[i,3]*H_XXYYZZ[:,:,3]
        # vals = eigvals(HamJij)
        # Heigenvals_f[:,i] = vals

        Sx_eq_rand_list = zeros(nrand_num)

        for j in 1:nrand_num
            rand_Jij = (randn() * sigmasqrt)
            Ham1 = rand_Jij * HamJ_noRandJ + Hamh0
            rho = exp(-beta*Ham1)
            rho = rho / tr(rho)

            # Sx per single site, and need 2 for 
            Sx_eq_rand_list[j] = real(tr(rho*Sx_op)) / N
        end
        Sx_eq_aver = mean(Sx_eq_rand_list)

        Sx_flist_val[i] = Sx_eq_aver
    end


    return Sx_flist_val
end


function main()
    
    # N = 6
    # disorder_num = 100

    N = 7
    nrand_num = 500
    # 
    # N = 6
    # disorder_num = 50
    # 
    # t_list only use the Delta_t


    data_dirname = "rs_Sx_phase_diagram_20220705"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end

    hinJ = 0.05
    J = 1.0


    # Deltalist = [-0.5, 0.0, 0.5]
    # betalist = [0.0, 4.0, 15.0, 20.0]
    # betalist = [4.0, 30.0]
    betalist = [0.1]
    # shorttime
    # t_list = collect(0:0.05:50)
    # longtime
    # t_list = collect(0:0.5:30)

    # Tlist = [10.0]

    # flist = [[1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]
    # xi_len must be even number
    # xi_len = 21
    xi_len = 2*10 + 1
    xi_step = 0.5

    flist = zeros(xi_len*xi_len,3)
    # ij_map = zeros(Int64, xi_len*xi_len, 2)
    # initial flist_len
    flist_len = 1
    for i in 1:xi_len
        for j in 1:xi_len
            flist[flist_len, 1] = 1
            flist[flist_len, 2] = (i - (xi_len - 1)÷2 -1) * xi_step
            flist[flist_len, 3] = (j - (xi_len - 1)÷2 -1) * xi_step
            
            # ij_map[flist_len, 1] = i
            # ij_map[flist_len, 2] = j

            flist_len += 1
        end
    end
    flist_len -= 1
    Jbar = 0.0

    Sx_flist_beta = zeros(flist_len, length(betalist))

    for (j, beta) in enumerate(betalist)
        # Delta = -0.73
        # Jbar = 0.0
        println("β=$(beta)")
        
        Sx_flist_val = random_realization_Sx(beta, N, J, Jbar,  hinJ, nrand_num, flist, flist_len)
        Sx_flist_beta[:,j] = Sx_flist_val
    end

    JLD.save("$(data_dirname)/Sx_flist_beta.jld", "Sx_flist_beta", Sx_flist_beta, "flist", flist)

    writedlm("$(data_dirname)/f_arr.csv", flist[:,:], ",")
    writedlm("$(data_dirname)/Sx_flist_beta.csv", Sx_flist_beta[:,:], ",")

    # JLD.save("$(data_dirname)/Sx_aver.jld", "Sx_aver", Sx_aver_listall, "t_list", tperbeta_list*betalist[1])
    # 
    # 
    # Sx_aver_listall = JLD.load("Sx_aver.jld", "Sx_aver")






    #################################################################################
    #=
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

    # JLD.save("$(data_dirname)/Sx_aver.jld", "Sx_aver", Sx_aver_listall)
    savefig("$(data_dirname)/thermalSx.pdf")
    =#

    ###################################### END Dynamics ###############################
end

main()