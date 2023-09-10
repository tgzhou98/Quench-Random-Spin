using Distributed
using JLD
using DelimitedFiles
using Plots
using LaTeXStrings
# gr()
# using Arpack
ENV["GKSwstype"]="100"
gr()

# # using ProgressMeter
# using Distributed
# using SharedArrays

include("rd_model_ed.jl")

function singleRand_mb_spectrum_func(N::Int64, J::Float64, Jbar::Float64, f::Array{Float64,1})
    # beta = 15
    # N = 8
    # println(f)
    HamJij = generate_RD_ham_Jij(N, J, Jbar, f)

    vals, vecs = eigen(HamJij)

    return vals
end


function singleRand_mb_spectrum_allf_func(N::Int64, J::Float64, Jbar::Float64, flist_arr::Array{Float64,2}, f_num::Int64)
    # beta = 15
    # N = 8
    # println(f)
    # HamJij = generate_RD_ham_Jij(N, J, Jbar, f)

    Heigenvals_f = zeros(2^N, f_num)
    
    H_XXYYZZ = generate_RD_ham_Jij_pauliXXYYZZ(N, J, Jbar)
    for i in 1:f_num
        println("N=$(N), J=$J, f=$(flist_arr[i,:]),  step=$i of $(f_num)")
        HamJij = flist_arr[i,1]*H_XXYYZZ[:,:,1] + flist_arr[i,2]*H_XXYYZZ[:,:,2] + flist_arr[i,3]*H_XXYYZZ[:,:,3]
        vals = eigvals(HamJij)
        Heigenvals_f[:,i] = vals
    end

    return Heigenvals_f
end

function random_realization_Heigenvals(Heigenvals_f::Array{Float64, 2}, nrand_num::Int64)
    sigmasqrt = sqrt(4*J^2/N)
    Hdim, f_num = size(Heigenvals_f)
    Heigenvals_f_rand = zeros(Hdim*nrand_num, f_num)
    for i in 1:nrand_num
        temp_rand = (randn() * sigmasqrt)
        Heigenvals_f_rand[1+(i-1)*Hdim:Hdim+(i-1)*Hdim, :] = temp_rand * Heigenvals_f
    end

    return Heigenvals_f_rand
end



N = 8
J = 1.
nrand_num = 200
# t_list only use the Delta_t


data_dirname = "rs_MB_spectrum_random_20220627"
if !isdir(data_dirname)
    mkdir(data_dirname)
end


# Deltalist = [0.0, -2.0, 2.0]
flist_in = [[1.0, -0.5, 0.0], [1.0, -0.8, -0.5], [1.0, -1.0, -1.0], [1.0, 1.0, 0.5], [1.0, 1.0, -0.5], [1., 1, -2], [1., 1, -1], [1., 1, 5], [1., 1, 10], [1., 1, -1.5], [1, 1, -2.5], [1, 1, -3.0], [1.0, -2.0, 0.0], [1.0, -1.2, 0.5], [1, 1.5, -0.5], [1.0, 0.8, -0.8], [1.0, -0.1, -2.0]]
f_num = length(flist_in)
flist = zeros(f_num, 3) 
for i in 1:f_num
    flist[i,:] = flist_in[i]
end

mb_spectrum_arr = zeros(Float64, 2^N, f_num)

Jbar = 0.0
#=
for i in 1:f_num
    println("N=$(N), J=$J, f=$(flist_in[i]),  step=$i of $(f_num)")
# for Jbar in [-0.3]
    f = flist[i,:]
    mb_spectrum_data = single_mb_spectrum_func(N, 1.0, Jbar, f)
    mb_spectrum_arr[:,i] = mb_spectrum_data

    JLD.save("$(data_dirname)/mb_spectrum.jld", "mb_spectrum_arr", mb_spectrum_arr)
end
=#
mb_spectrum_arr_single = singleRand_mb_spectrum_allf_func(N, J, Jbar, flist, f_num)
mb_spectrum_arr = random_realization_Heigenvals(mb_spectrum_arr_single, nrand_num)


writedlm("$(data_dirname)/f_arr.csv", flist[:,:], ",")
writedlm("$(data_dirname)/mb_spectrum_arr.csv", mb_spectrum_arr, ",")
# writedlm("$(data_dirname)/t_list.csv", tperbeta_list*betalist[1], ",")

# JLD.save("$(data_dirname)/Sx_aver.jld", "Sx_aver", Sx_aver_listall, "t_list", tperbeta_list*betalist[1])
# 
# 
# Sx_aver_listall = JLD.load("Sx_aver.jld", "Sx_aver")


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

# JLD.save("$(data_dirname)/Sx_aver.jld", "Sx_aver", Sx_aver_listall)
savefig("$(data_dirname)/thermalSx.pdf")
# savefig("thermalSx.pdf")

# Sx_aver_listall = JLD.load("$(data_dirname)/Sx_aver.jld", "Sx_aver")
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
for ll in 1:f_num
    label_f = "\$ f=$(flist[ll,:]) \$"
    A = -(flist[ll,1])^2 + (flist[ll,2])^2 - 4*flist[ll,2]*flist[ll,3] + (flist[ll,3])^2

    # p1 = histogram(mb_spectrum_arr[:,ll], title = "\$ \\beta=$(betalist[1]), N=$(N), \\bar{J}=$(Jbar), h_i=$(hinJ), h_f=$(quench_hinJ) \$", xlabel="\$ tJ \$",ylabel=L" \textrm{Tr } (\rho \hat{S}^x(t)) ", label=label_f, titlefont=9, legendfontsize=7, xguidefontsize=7,yguidefontsize=7)
    histogram(mb_spectrum_arr[:,ll], title = "\$ N=$(N), J=$(J), \\bar{J}=$(Jbar), f=$(flist[ll,:]), A(f)=$(A)\$", nbins=300, xlabel="\$ E \$",ylabel=L"$N(E)$")

    # JLD.save("$(data_dirname)/Sx_aver.jld", "Sx_aver", Sx_aver_listall)
    savefig("$(data_dirname)/manybody_spectrum_$(ll).png")
    savefig("$(data_dirname)/manybody_spectrum_$(ll).pdf")
end



###################################### END Dynamics ###############################
