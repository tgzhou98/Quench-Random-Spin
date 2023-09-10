include("miscellaneous.jl")
# Warning Also load miscellaneous.jl first because it will be use on other file. Not to reinclude in other file!!!
include("equilibrium.jl")
include("SY_dynamics.jl")
using Glob

function get_dynmaics(βJ::Float64, hinJ::Float64, field::String, Δ::Float64, Λω::Float64, Λt::Float64, nt::Int64, init_file::String, mode::String)
    Δt = Λt/((nt-1)÷2)


    ################################################################
    # equilibrium
    ################################################################
    ρω, abscissaeω=equilibrium(βJ, hinJ, field, Δ, Λω, Λt, init_file)


    ################################################################
    # dynamics
    ################################################################

    hinJ_real = 0.0
    if mode == "TI"
        hinJ_real = hinJ
    elseif mode == "quench"
        hinJ_real = 0.0
    end
    J, H0, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
     dynamics_init(βJ, hinJ_real, field, Δ, ρω, Λω, Λt, nt);

    baym_kadanoff(J, H0, nt, Δt, f
        ,GGreaterDynamicsttp  
        ,GLesserDynamicsttp   
        ,ΣGreaterDynamicsttp  
        ,ΣLesserDynamicsttp   
        ,GRetardedDynamicsttp 
        ,GAdavancedDynamicsttp
        ,ΣRetardedDynamicsttp 
        ,ΣAdavancedDynamicsttp
        ,mode
        ,false
    )

    # temp_file = "./quench_data/betaJ9_0/Dynamics.jld"
    # GGreaterDynamicsttp = JLD.load(temp_file, "GGreaterDynamicsttp");
    # GLesserDynamicsttp= JLD.load(temp_file, "GLesserDynamicsttp");

    # Plot after quench
    center_pt = (nt+1)÷2
    plot_Gt(GLesserDynamicsttp,   βJ, hinJ, field, Δ, β, center_pt - 200, center_pt + (nt-1)÷2÷4, Δt, center_pt)
    plot_Sx(GLesserDynamicsttp,   βJ, hinJ, field, Δ, β, center_pt - 200, nt, Δt, center_pt)
    plot_n_up(GLesserDynamicsttp, βJ, hinJ, field, Δ, β, center_pt - 200, nt, Δt, center_pt)
    export_Sx(GLesserDynamicsttp)

    # Move file 
    data_dirname = "quench_data"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    betaJ_string = replace("$(βJ)","."=>"_")
    target_dir = "$(data_dirname)/betaJ$(betaJ_string)"
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    for jld_file in ["Dynamics.jld", "ρω.jld"]
        if isfile(jld_file)
            mv(jld_file, "$(target_dir)/$(jld_file)", force=true)
        end
    end
    for png_file in glob("*.png")
        mv(png_file, "$(target_dir)/$(png_file)", force=true)
    end
    for pdf_file in glob("*.pdf")
        mv(pdf_file, "$(target_dir)/$(pdf_file)", force=true)
    end
    for csv_file in ["Sx.csv"]
        if isfile(csv_file)
            mv(csv_file, "$(target_dir)/$(csv_file)", force=true)
        end
    end
end

function get_ρω(βJ::Float64, hinJ::Float64, field::String, Δ::Float64, Λω::Float64, Λt::Float64, init_file::String)
    ################################################################
    # equilibrium
    ################################################################
    ρω, abscissaeω=equilibrium(βJ, hinJ, field, Δ, Λω, Λt, init_file)

    # Move plot and jld file saving to equilibrium_Jbar.jl

    # Move file 
    data_dirname = "spectrum_anal"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    betaJ_string = replace("$(βJ)","."=>"_")
    target_dir = "$(data_dirname)/betaJ$(betaJ_string)"
    if !isdir(target_dir)
        mkdir(target_dir)
    end

    for ρω_file in glob("ρω*.*")
        mv(ρω_file, "$(target_dir)/$(ρω_file)", force=true)
    end
    return ρω, abscissaeω

    # # Transfer to MATLAB
    # data_dirname = "spectrum_anal"
    # if !isdir(data_dirname)
    #     mkdir(data_dirname)
    # end
    # betaJ_string = replace("$(βJ)","."=>"_")
    # target_dir = "$(data_dirname)/betaJ$(betaJ_string)"
    # ρω = JLD.load("$(target_dir)/ρω0.jld", "arr")
    # file = matopen("$(target_dir)/rhoomega.mat", "w")
    # write(file, "rhoomega", ρω)
    # close(file)

    # file = matopen("omega.mat", "w")
    # write(file, "omega", abscissaeω)
    # close(file)

    # return ρω
end

function run_equilibrium(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, hinJ::Float64, field::String, Δ::Float64, Λω::Float64, Λt::Float64, init_file::String)

    run_num = Int64((end_betaJ - begin_betaJ)/step_betaJ) + 1
    ρωlist = Array{ComplexF64, 4}(undef, run_num, 3003, 2, 2)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    abscissaeω = Array{ComplexF64, 3}(undef, 0, 0, 0)

    # for i in 4:14
    for i in 1:run_num
        ρωlist[i, :,:,:], abscissaeω = get_ρω(betaJlist[i], hinJ, field, Δ, Λω, Λt, init_file)
    end

    # TODO, analyze the spectrum information
    # peak_analyze(run_num, ρωlist, abscissaeω, betaJlist)
end



function run_dynmaics(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, hinJ::Float64, field::String, Δ::Float64, Λω::Float64, Λt::Float64, nt::Int64, init_file::String, mode::String = "quench", fix_h::Bool = false, fix_h_betaJ::Float64 = 6.0)

    run_num = Int64((end_betaJ - begin_betaJ)/step_betaJ) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)

    if fix_h
        for i in 1:run_num
            get_dynmaics(betaJlist[i], hinJ*fix_h_betaJ/betaJlist[i], field, Δ, Λω, Λt, nt, init_file, mode)
        end
    else
        for i in 1:run_num
            get_dynmaics(betaJlist[i], hinJ, field, Δ, Λω, Λt, nt, init_file, mode)
        end
    end
end

function peaks_analyze()
    Λω = 10

    abscissae1 = JLD.load("IntegrationData.jld", "abscissae1")
    weights1= JLD.load("IntegrationData.jld", "weights1")
    errweights1= JLD.load("IntegrationData.jld", "errweights1")

    nω = length(abscissae1)
    a1 = -Λω + 0.00000001
    b1 = Λω
    abscissaeω = [(b1 - a1)*abscissae1[i] + a1 for i in 1:nω]
    weightsω = [(b1 - a1)*weights1[i] for i in 1:nω]

    abscissaeω = vcat((abscissaeω[(1 + nω)÷2:nω] .- Λω), [0.00000000001], (Λω .+ abscissaeω[1 : (1 + nω)÷2]))
    weightsω =   vcat((weightsω[(1 + nω)÷2:nω]), [0.], (weightsω[1 : (1 + nω)÷2]))

    for betaJ in 16:21
        ρω=JLD.load("./quench/betaJ$(betaJ)_0/ρω.jld","arr")
        peaks_pt = findlocalmaxima(real(ρω[:,1,1]))
        peaks_ω = [abscissaeω[i] for i in peaks_pt] 
        # println(peaks)
        # println([real(ρω[i,1,1]) for i in peaks])
        @assert length(peaks)==4
        nω = size(ρω)[1]
        center_pt = (nω+1)÷2
        ralative_peaks = abs.(peaks_ω .- abscissaeω[center_pt])
        ratio_peaks12 = (ralative_peaks[1]/ralative_peaks[2] + ralative_peaks[4]/ralative_peaks[3])/2
        
        println("βJ=$(betaJ), peaks ratio = $(ratio_peaks12)")
    end
end

# peaks_analyze()

######################################################################################################
# h not fixed
######################################################################################################
# # init_file = "./spectrum_anal_break_phase/betaJ12_0/ρω.jld"
# init_file = ""
# # run_equilibrium(13.0, 13.0, 1.0, 0.0, "x", -0.73, 15.0, 50.0, init_file)
# # run_dynmaics(2.0, 5.0, 1.0, 0.05, "x", -0.73, 10.0, 75.0, 4001, init_file, "TI")
# data_dirname = "quench_data_Delta_M0_73_high_accu"
# if !isdir(data_dirname)
#     mkdir(data_dirname)
# end
# run_dynmaics(2.0, 8.0, 1.0, 0.05, "x", -0.73, 8.0, 75.0, 4001, init_file, "quench")
# mv("quench_data", "$(data_dirname)/quench")
# run_dynmaics(2.0, 8.0, 1.0, 0.05, "x", -0.73, 8.0, 75.0, 4001, init_file, "TI")
# mv("quench_data", "$(data_dirname)/TI")

######################################################################################################
# h fixed
######################################################################################################
# init_file = "./spectrum_anal_break_phase/betaJ12_0/ρω.jld"
init_file = ""
# run_equilibrium(13.0, 13.0, 1.0, 0.0, "x", -0.73, 15.0, 50.0, init_file)
# run_dynmaics(2.0, 5.0, 1.0, 0.05, "x", -0.73, 10.0, 75.0, 4001, init_file, "TI")
data_dirname = "quench_data_Delta_M0_73_high_accu_fixh"
if !isdir(data_dirname)
    mkdir(data_dirname)
end
run_dynmaics(1.0, 6.0, 1.0, 1.25, "x", -0.73, 8.0, 75.0, 4001, init_file, "quench", true, 6.0)
mv("quench_data", "$(data_dirname)/quench")
run_dynmaics(1.0, 6.0, 1.0, 1.25, "x", -0.73, 8.0, 75.0, 4001, init_file, "TI", true, 6.0)
mv("quench_data", "$(data_dirname)/TI")

# init_file = ""
# run_dynmaics(2.0, 2.0, 1.0, 0.05, "x", -0.73, 8.0, 75.0, 2001, init_file, "TI")


######################################################################################################
# Single shot run in global scope
######################################################################################################
# ρω, abscissaeω=equilibrium(2.0, 0.05, "x", -0.73, 8.0, 50.0, "")
# J, H0, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
#     dynamics_init(2.0, 0.05, "x", -0.73, ρω, 8.0, 50.0, 2001);
# baym_kadanoff(J, H0, nt, Δt, f
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