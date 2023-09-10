include("miscellaneous.jl")
# Warning Also load miscellaneous.jl first because it will be use on other file. Not to reinclude in other file!!!
include("equilibrium_Jbar_NMR_H1_H2.jl")
include("SY_dynamics_Jbar.jl")
using Glob

function write_log(data_dirname::String, begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64 = 0.0, Λt_arr::Array{Float64, 2} = zeros(1,1))
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    open("$(data_dirname)/parameter.csv", "w") do io
        write(io,
        """
        begin_betaJ, $(begin_betaJ)
        end_betaJ, $(end_betaJ)
        step_betaJ, $(end_betaJ)
        begin_hinJ, $(begin_hinJ)
        end_hinJ, $(end_hinJ)
        step_hinJ, $(step_hinJ)
        std_betaJ, $(std_betaJ)
        field, $(field)
        f, $(f)
        JbarinJ, $(JbarinJ)
        Λω, $(Λω)
        Λt_arr, $(Λt_arr)
        """
        )
    end;
end


function write_log(data_dirname::String, Tlist::Array{Float64, 1},
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64 = 0.0, Λt_arr::Array{Float64, 2} = zeros(1,1))
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    open("$(data_dirname)/parameter.csv", "w") do io
        write(io,
        """
        Tlist, $(Tlist)
        begin_hinJ, $(begin_hinJ)
        end_hinJ, $(end_hinJ)
        step_hinJ, $(step_hinJ)
        std_betaJ, $(std_betaJ)
        field, $(field)
        f, $(f)
        JbarinJ, $(JbarinJ)
        Λω, $(Λω)
        Λt_arr, $(Λt_arr)
        """
        )
    end;
end

function phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff_generate(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
     begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt_arr::Array{Float64, 2}, nt::Int64,
     init_GRω::Array{ComplexF64, 3},
     mode::String = "quench", fix_beta::Bool = true)

    new_data_dirname = "phase_diagram_dynamics"
    write_log(new_data_dirname, begin_betaJ, end_betaJ, step_betaJ, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ, Λω, Λt_arr)
    export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)
end

function run_dynmaics_Jbar(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, hinJ::Float64, field::String, f::Array{Float64,1}, JbarinJ::Float64,
        Λω::Float64, Λt::Float64, nt::Int64, mode::String = "quench", fix_h::Bool = false, fix_h_betaJ::Float64 = 6.0, fix_beta::Bool = true, J_std::Float64 = 1.0)

    run_num = Int64((end_betaJ - begin_betaJ)/step_betaJ) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)

    # TODO I don't want to write this right now
    β = 2π

    if fix_h
        for i in 1:run_num
            get_dynmaics_Jbar(betaJlist[i], hinJ*fix_h_betaJ/betaJlist[i], field, f, JbarinJ, β, Λω, Λt, nt, mode)
        end
    else
        for i in 1:run_num
            get_dynmaics_Jbar(betaJlist[i], hinJ, field, f, JbarinJ, β, Λω, Λt, nt, mode)
        end
    end
end

function get_dynmaics_Jbar(βJ::Float64, hinJ::Float64, field::String, f::Array{Float64,1}, JbarinJ::Float64, β::Float64,
    Λω::Float64, Λt::Float64, nt::Int64, mode::String, init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0))
    Δt = Λt/((nt-1)÷2)


    ################################################################
    # equilibrium
    ################################################################

    println("---------------------------------")
    println("equilibrium SCF Start: T=$(1/β), betaJ=$(βJ), h/J=$(hinJ), f=$(f), Λt=$(Λt), Λω=$(Λω)")
    iter_mode = false
    ρω, abscissaeω, Sxyz_new, h0, Jbar, GRω = equilibrium_Jbar(βJ, hinJ, field, f, JbarinJ, β, Λω , Λt, init_GRω, iter_mode)

    spectral_ratio_val = spectral_ratio(ρω)
    println("Spectral Ratio=$(spectral_ratio_val)")
    # beta = 2π
    println("h0=$(h0)")
    println("h=$(hinJ*βJ/(β))", ",  Jbar=$(Jbar)")
    println("Sx:$(Sxyz_new[1]), Sy:$(Sxyz_new[2]), Sz:$(Sxyz_new[3]), sqrt(S^2):$(Sxyz_new[4])")
        # dirname_temp = "quench_data_Jbar_Delta_M0_73_high_accu_thermal_rhoomega/Jbar1_5_cutoff_75/betaJ0_5/"
    # dirname_temp = "./quench_data_Jbar/betaJ0_5/"
    # dirname_temp = "./phase_diagram_dynamics/betaJ0_5/"
    
    # ρω = JLD.load("$(dirname_temp)ρω.jld", "arr")
    # abscissaeω = JLD.load("$(dirname_temp)ρω.jld", "omega")
    # Sx_new = JLD.load("$(dirname_temp)ρω.jld", "Sx_new")
    # h0 = JLD.load("$(dirname_temp)ρω.jld", "h0")
    # Jbar = JLD.load("$(dirname_temp)ρω.jld", "Jbar")

    println("h0=$(h0)")
    println("Jbar=$(Jbar)")
    println("Sx=$(Sxyz_new)")
    println("βJ=$(βJ)")

    ################################################################
    # dynamics
    ################################################################

    J, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
        dynamics_init_Jbar(βJ, field, f, β, ρω, Λω, Λt, nt);

    showstep = 1
    println("###################################################")
    println("Start printing initial solution")
    println("G_>(t1,t2)=")
    display(GGreaterDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("G_<(t1,t2)=")
    display(GLesserDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("G_K(t1,t2)=")
    display((GLesserDynamicsttp+GGreaterDynamicsttp)[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_>(t1,t2)=")
    display(ΣGreaterDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_<(t1,t2)=")
    display(ΣLesserDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_K(t1,t2)=")
    display((ΣLesserDynamicsttp+ΣGreaterDynamicsttp)[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("###################################################")

    hinJ_real = 0.0
    if mode == "TI"
        hinJ_real = hinJ
    elseif mode == "quench"
        hinJ_real = 0.0
    end

    k_index = baym_kadanoff_Jbar(J, hinJ_real, field, Sxyz_new, Jbar, nt, Δt, f
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
    println("end index: ", k_index)
    ΣKDynamicsttp = ΣGreaterDynamicsttp + ΣLesserDynamicsttp

    println("###################################################")
    println("Start printing converged solution")
    println("G_>(t1,t2)=")
    display(GGreaterDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("G_<(t1,t2)=")
    display(GLesserDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("G_K(t1,t2)=")
    display((GLesserDynamicsttp+GGreaterDynamicsttp)[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_>(t1,t2)=")
    display(ΣGreaterDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_<(t1,t2)=")
    display(ΣLesserDynamicsttp[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("Σ_K(t1,t2)=")
    display((ΣLesserDynamicsttp+ΣGreaterDynamicsttp)[(nt+1)÷2-showstep:(nt+1)÷2+showstep,(nt+1)÷2-showstep:(nt+1)÷2+showstep,:,:])
    println("###################################################")


    # Plot after quench
    center_pt = (nt+1)÷2
    plot_Gt(GLesserDynamicsttp,   βJ, hinJ, field, f, β, center_pt - 200, k_index, Δt, center_pt)
    plot_allS(GLesserDynamicsttp,   βJ, hinJ, JbarinJ, field, f, β, center_pt - 200, nt, Δt, center_pt)
    plot_n_up(GLesserDynamicsttp, βJ, hinJ, field, f, β, center_pt - 200, nt, Δt, center_pt)

    rho_list_num = 16
    # plot_deltak_index = nt*3÷10
    plot_deltak_index = k_index - center_pt - 5
    # ρωthermal_list = Array{ComplexF64, 4}(undef, length(abscissaeω), 2, 2, rho_list_num)
    plot_t_list = Array{Float64, 1}(undef, rho_list_num + 1)
    cutoff_length = (nt - 1) ÷ 2

    # get t list
    # Λt = 75
    a2=-Λt+0.00000001;
    b2= Λt;
    # nt = 2001;
    equaldTabscissae = collect(0:1/(nt-1):1)
    Δt=Λt/((nt-1)/2);
    equaldTabscissaet=[(b2-a2)*equaldTabscissae[i]+a2 for i in 1:nt];


    for i in 1:rho_list_num+1
        plot_thermal_step = plot_deltak_index*(i-1)÷(rho_list_num) + center_pt
        plot_t_list[i] = (plot_thermal_step-center_pt)*Δt
        
        # Plot thermal spectrum
        ρωthermal_temp, Fomegathermal = thermal_spectrum_F(ΣKDynamicsttp, GRetardedDynamicsttp, ΣRetardedDynamicsttp, β, Λω, Λt, nt, plot_thermal_step)
        plot_ρω(ρωthermal_temp, abscissaeω, βJ, hinJ, field, f, β, (plot_thermal_step-center_pt)*Δt, i, "thermal")

        # Plot thermal Gt
        Gtthermal_temp, Glessomegathermal = thermal_Gt(GLesserDynamicsttp, β, Λω, Λt, nt, plot_thermal_step)
        plot_Gt_eq(Gtthermal_temp, equaldTabscissaet, βJ, hinJ, field, f, β, (plot_thermal_step-center_pt)*Δt, i, "thermal")
        save_arrayω("distF", Fomegathermal, abscissaeω, "real", i, "thermal")

    end
    
    writedlm("abscissaet_beta.csv", equaldTabscissaet/β, ',')
    writedlm("abscissaeω_J.csv", abscissaeω/J, ',')
    writedlm("plot_t_list.csv", plot_t_list, ',')

    # JLD.save("ρω_thermal.jld","ρωthermal_list", ρωthermal_list, "omega", abscissaeω, "plot_t_list", plot_t_list)
    export_Sx(GLesserDynamicsttp)

    # Move file 
    data_dirname = "quench_data_Jbar"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    # betaJ_string = replace("$(round(βJ, digits=2))","."=>"_")
    ToverJ_string = replace("$(round(1/βJ, digits=3))","."=>"_")
    hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
    # target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"
    target_dir = "$(data_dirname)/ToverJ$(ToverJ_string)_hinJ$(hinJ_string)"
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    for jld_file in ["Dynamics.jld", "ρω.jld", "ρω_thermal.jld"]
        if isfile(jld_file)
            mv(jld_file, "$(target_dir)/$(jld_file)", force=true)
        end
    end
    for pdf_file in glob("*.pdf")
        mv(pdf_file, "$(target_dir)/$(pdf_file)", force=true)
    end
    for png_file in glob("*.png")
        mv(png_file, "$(target_dir)/$(png_file)", force=true)
    end
    for csv_file in glob("*.csv")
        if isfile(csv_file) && csv_file != "IntegrationData.csv"
            mv(csv_file, "$(target_dir)/$(csv_file)", force=true)
        end
    end
    mv("GLessert.jld", "$(target_dir)/GLessert.jld", force=true)

    # Sxlist_t = imag(diag(GLesserDynamicsttp[:,:,1,2])[center_pt:nt])

    return GRω
    
end

# initial x: hinJ
# evolution z: hinJ/4
function get_dynmaics_Jbar_fieldxtoz(βJ::Float64, hinJ::Float64, f::Array{Float64,1}, JbarinJ::Float64, β::Float64,
    Λω::Float64, Λt::Float64, nt::Int64, mode::String, init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0))
    Δt = Λt/((nt-1)÷2)


    ################################################################
    # equilibrium
    ################################################################
    initial_field = "x"
    ρω, abscissaeω, Sxyz_new, h0, Jbar, GRω = equilibrium_Jbar(βJ, hinJ, initial_field, f, JbarinJ, β, Λω , Λt, init_GRω)
        # dirname_temp = "quench_data_Jbar_Delta_M0_73_high_accu_thermal_rhoomega/Jbar1_5_cutoff_75/betaJ0_5/"
    # dirname_temp = "./quench_data_Jbar/betaJ0_5/"
    # dirname_temp = "./phase_diagram_dynamics/betaJ0_5/"
    
    # ρω = JLD.load("$(dirname_temp)ρω.jld", "arr")
    # abscissaeω = JLD.load("$(dirname_temp)ρω.jld", "omega")
    # Sx_new = JLD.load("$(dirname_temp)ρω.jld", "Sx_new")
    # h0 = JLD.load("$(dirname_temp)ρω.jld", "h0")
    # Jbar = JLD.load("$(dirname_temp)ρω.jld", "Jbar")

    println("h0=$(h0)")
    println("Jbar=$(Jbar)")
    println("Sx:$(Sxyz_new[1]), Sy:$(Sxyz_new[2]), Sz:$(Sxyz_new[3]), sqrt(S^2):$(Sxyz_new[4])")
    println("x to z experiment mode")

    ################################################################
    # dynamics
    ################################################################

    quench_field = "z"
    J, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp =
        dynamics_init_Jbar(βJ, quench_field, f, β, ρω, Λω, Λt, nt);

    # hinJ_real = 0.0
    # if mode == "TI"
    #     hinJ_real = hinJ
    # elseif mode == "quench"
    #     hinJ_real = 0.0
    # end

    k_index = baym_kadanoff_Jbar(J, hinJ/8, quench_field, Sxyz_new, Jbar, nt, Δt, f
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
    println("end index: ", k_index)

    # Plot after quench
    center_pt = (nt+1)÷2
    plot_Gt(GLesserDynamicsttp,   βJ, hinJ, quench_field, f, β, center_pt - 200, k_index, Δt, center_pt)
    # plot_Sx(GLesserDynamicsttp,   βJ, hinJ, quench_field, f, β, center_pt - 200, nt, Δt, center_pt)
    plot_allS(GLesserDynamicsttp,   βJ, hinJ, JbarinJ, quench_field, f, β, center_pt - 200, nt, Δt, center_pt)
    plot_n_up(GLesserDynamicsttp, βJ, hinJ, quench_field, f, β, center_pt - 200, nt, Δt, center_pt)

    rho_list_num = 20
    plot_deltak_index = 800
    ρωthermal_list = Array{ComplexF64, 4}(undef, length(abscissaeω), 2, 2, rho_list_num)
    plot_t_list = Array{Float64, 1}(undef, rho_list_num)
    for i in 1:rho_list_num
        plot_thermal_step = plot_deltak_index*i÷rho_list_num + center_pt
        plot_t_list[i] = (plot_thermal_step-center_pt)*Δt
        ρωthermal_list[:,:,:,i] = thermal_spectrum_F(GRetardedDynamicsttp, β, Λω, Λt, nt, plot_thermal_step)
        plot_ρω(ρωthermal_list[:,:,:,i], abscissaeω, βJ, hinJ, quench_field, f, β, (plot_thermal_step-center_pt)*Δt, i, "thermal")
    end
    JLD.save("ρω_thermal.jld","ρωthermal_list", ρωthermal_list, "omega", abscissaeω, "plot_t_list", plot_t_list)
    export_Sx(GLesserDynamicsttp)

    # Move file 
    data_dirname = "quench_data_Jbar"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    betaJ_string = replace("$(round(βJ, digits=2))","."=>"_")
    hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
    target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    for jld_file in ["Dynamics.jld", "ρω.jld", "ρω_thermal.jld"]
        if isfile(jld_file)
            mv(jld_file, "$(target_dir)/$(jld_file)", force=true)
        end
    end
    for pdf_file in glob("*.pdf")
        mv(pdf_file, "$(target_dir)/$(pdf_file)", force=true)
    end
    for png_file in glob("*.png")
        mv(png_file, "$(target_dir)/$(png_file)", force=true)
    end
    for csv_file in ["Sx.csv"]
        if isfile(csv_file)
            mv(csv_file, "$(target_dir)/$(csv_file)", force=true)
        end
    end
    mv("GLessert.jld", "$(target_dir)/GLessert.jld", force=true)

    # Sxlist_t = imag(diag(GLesserDynamicsttp[:,:,1,2])[center_pt:nt])

    return GRω
    
end


function get_ρω_Jbar(βJ::Float64, hinJ::Float64, field::String, f::Array{Float64,1}, JbarinJ::Float64, β::Float64, Λω::Float64, Λt::Float64, init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0), iter_mode::Bool = false, TorInvT::Bool = true)
    ################################################################
    # equilibrium
    ################################################################
    println("---------------------------------")
    println("equilibrium SCF Start: T=$(1/β), betaJ=$(βJ), h/J=$(hinJ), f=$(f), Jbar/J=$(JbarinJ)")
    ρω, abscissaeω, Sxyz_new, h0, Jbar, GRω = equilibrium_Jbar(βJ, hinJ, field, f, JbarinJ, β, Λω , Λt, init_GRω, iter_mode)
    # beta = 2π
    println("h0=$(h0)")
    println("h=$(hinJ*βJ/(β))", ",  Jbar=$(Jbar)")
    println("Sx:$(Sxyz_new[1]), Sy:$(Sxyz_new[2]), Sz:$(Sxyz_new[3]), sqrt(S^2):$(Sxyz_new[4])")

    # Move plot and jld file saving to equilibrium_Jbar.jl

    # Move file 
    data_dirname = "spectrum_anal_Jbar"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    # betaJ_string = replace("$(round(βJ, digits=2))","."=>"_")
    # hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
    # target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"

    if TorInvT
        T_string = replace("$(round(1/β, digits=6))","."=>"_")
        hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
        target_dir = "$(data_dirname)/T$(T_string)_hinJ$(hinJ_string)"
    else
        betaJ_string = replace("$(round(βJ, digits=6))","."=>"_")
        hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
        target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"
    end

    if !isdir(target_dir)
        mkdir(target_dir)
    end
    for ρω_file in glob("ρω*.*")
        mv(ρω_file, "$(target_dir)/$(ρω_file)", force=true)
    end
    for Gt_file in glob("Gt*.*")
        mv(Gt_file, "$(target_dir)/$(Gt_file)", force=true)
    end
    for GRω_file in glob("GRω*.*")
        mv(GRω_file, "$(target_dir)/$(GRω_file)", force=true)
    end
    for pdf_file in glob("*.pdf")
        mv(pdf_file, "$(target_dir)/$(pdf_file)", force=true)
    end
    for png_file in glob("*.png")
        mv(png_file, "$(target_dir)/$(png_file)", force=true)
    end
    mv("GLessert.jld", "$(target_dir)/GLessert.jld", force=true)
    
    return ρω, Sxyz_new, h0, GRω
end


function get_Gt_imagtime_Jbar(βJ::Float64, hinJ::Float64, field::String, f::Array{Float64,1}, JbarinJ::Float64, β::Float64,
    init_Gtau::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0), TorInvT::Bool = false, SelfEtype::String = "std", iter_mode::Bool = false)
    ################################################################
    # equilibrium
    ################################################################
    println("---------------------------------")
    println("equilibrium SCF Start: T=$(1/β), Jbar/J=$(JbarinJ), betaJ=$(βJ), h/J=$(hinJ), f=$(f)")

    abscissaet, Sxyz_new, h0, free_energy, Jbar, Gt = equilibrium_Jbar_imagtime(βJ, hinJ, field, f, JbarinJ, β, init_Gtau, SelfEtype, iter_mode)
    # beta = 2π
    println("h0=$(h0)")
    println("h=$(hinJ*βJ/(β))", ",  Jbar=$(Jbar)")
    println("Sx:$(Sxyz_new[1]), Sy:$(Sxyz_new[2]), Sz:$(Sxyz_new[3]), sqrt(S^2):$(Sxyz_new[4])")

    # Move plot and jld file saving to equilibrium_Jbar.jl

    # Move file 
    data_dirname = "Imagtime_eq_anal_Jbar"
    if !isdir(data_dirname)
        mkdir(data_dirname)
    end
    if TorInvT
        T_string = replace("$(round(1/β, digits=6))","."=>"_")
        hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
        target_dir = "$(data_dirname)/T$(T_string)_hinJ$(hinJ_string)"
    else
        betaJ_string = replace("$(round(βJ, digits=6))","."=>"_")
        hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
        target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"
    end
    if !isdir(target_dir)
        mkdir(target_dir)
    end
    for imagtime_eq_file in glob("imagtime_eq*.*")
        mv(imagtime_eq_file, "$(target_dir)/$(imagtime_eq_file)", force=true)
    end
    for csv_file in glob("*.csv")
        if isfile(csv_file) && csv_file != "IntegrationData.csv"
            mv(csv_file, "$(target_dir)/$(csv_file)", force=true)
        end
    end
    for pdf_file in glob("*.pdf")
        mv(pdf_file, "$(target_dir)/$(pdf_file)", force=true)
    end
    for png_file in glob("*.png")
        mv(png_file, "$(target_dir)/$(png_file)", force=true)
    end
    # mv("GLessert.jld", "$(target_dir)/GLessert.jld", force=true)
    
    return Sxyz_new, h0, free_energy, Gt
end


function run_equilibrium_Jbar(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, hinJ::Float64, field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt::Float64, β::Float64 = 2π,  init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0))

    run_num = Int64((end_betaJ - begin_betaJ)/step_betaJ) + 1
    ρωlist = Array{ComplexF64, 4}(undef, run_num, 3003, 2, 2)
    Sxyzlist = Array{Float64, 2}(undef, run_num, 4)
    h0list = Array{Float64, 1}(undef, run_num)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)

    # for i in 4:14
    for i in 1:run_num
        ρωlist[i, :,:,:], Sxyzlist[i, :], h0list[i], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJ, field, f, β, Λω, Λt)
    end

    # TODO, analyze the spectrum information
    return ρωlist, Sxyzlist, h0list
    
end


function run_phase_diagram_dynamics_Jbar_fixβ_J(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt::Float64, nt::Int64, mode::String = "quench", fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    # Sx_t_list = Array{Float64, 3}(undef, (nt + 1)÷2, betaJ_num, h_num)

    if fix_beta
        println("\n\nfix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for i in 1:betaJ_num
            for j in 1:h_num
               get_dynmaics_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, nt, mode)
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    else
        # println("fix beta")
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for i in 1:betaJ_num
            for j in 1:h_num
                get_dynmaics_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, nt, mode)
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist, "Sx_t_list")
    end

    old_data_dirname = "quench_data_Jbar"
    new_data_dirname = "phase_diagram_dynamics"
    if isdir(new_data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram_dynamics.jld", "$(new_data_dirname)/phase_diagram_dynamics.jld")
        println("file moved")
    end
end


function run_phase_diagram_Jbar_fixβ_J(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt::Float64, fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    ρω_centerlist = Array{ComplexF64, 4}(undef, betaJ_num, h_num, 2, 2)
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    if fix_beta
        println("fix beta")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for i in 1:betaJ_num
            for j in 1:h_num
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j] = get_ρω_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt)
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for i in 1:betaJ_num
            for j in 1:h_num
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j] = get_ρω_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt)
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    end

    old_data_dirname = "spectrum_anal_Jbar"
    new_data_dirname = "phase_diagram"
    if isdir(new_data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram.jld", "$(new_data_dirname)/phase_diagram.jld")
        println("file moved")
    end

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end

function run_plot_phase_diagram(dirname::String, filename::String, fix_beta::Bool = true)
    path_name = "$(dirname)/$(filename)"

    betaJlist = JLD.load(path_name, "betaJlist")
    hlist = JLD.load(path_name, "hlist")
    ρω_centerlist = JLD.load(path_name, "ρω_centerlist")
    Sxyzlist = JLD.load(path_name, "Sxyzlist")
    # Sxlist = JLD.load(path_name, "Sxlist")
    h0list = JLD.load(path_name, "h0list")
    # println(real(ρω_centerlist[:,:,1,1]))
    # println(round.(hlist*2π, digits=2))

    plot_phase_diagram(betaJlist, hlist, real(ρω_centerlist[:,:,1,1]), Sxyzlist[:,1], h0list, fix_beta)

    for png_file in glob("phase_diagram*.pdf")
        mv(png_file, "$(dirname)/$(png_file)", force=true)
    end
end

# Note that Λω::Float64 doesn't have to scale
# this is modified from Λω_scale mode. I don't need this
function run_phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
     begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, scale_Λt::Float64, nt::Int64,
     init_GRω::Array{ComplexF64, 3},
     mode::String = "quench", fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    # Sx_t_list = Array{Float64, 3}(undef, (nt + 1)÷2, betaJ_num, h_num)

    nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    if fix_beta
        println("fix beta")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    end

    old_data_dirname = "quench_data_Jbar"
    new_data_dirname = "phase_diagram_dynamics"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram_dynamics.jld", "$(new_data_dirname)/phase_diagram_dynamics.jld")
        println("file moved")
    end
end

# this is modified from Λω_scale mode. I don't need this
function run_phase_diagram_Jbar_fixβ_J_tune_init_GRω(begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64, begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, scale_Λt::Float64, init_GRω::Array{ComplexF64, 3}, fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    ρω_centerlist = Array{ComplexF64, 4}(undef, betaJ_num, h_num, 2, 2)
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    if fix_beta
        println("fix beta")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]

                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]

                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    end

    old_data_dirname = "spectrum_anal_Jbar"
    new_data_dirname = "phase_diagram"
    if isdir(data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram.jld", "$(new_data_dirname)/phase_diagram.jld")
        println("file moved")
    end

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end

function run_phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
     begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt_arr::Array{Float64, 2}, nt::Int64,
     init_GRω::Array{ComplexF64, 3},
     mode::String = "quench", fix_beta::Bool = true, nω::Int64 = 3003)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    # Sx_t_list = Array{Float64, 3}(undef, (nt + 1)÷2, betaJ_num, h_num)

    # nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    if fix_beta
        println("fix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                #####################################
                # Warning: fix_beta part may not work
                #####################################
                Λt = Λt_arr[i, j]
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    else
        # println("fix beta")
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Λt = Λt_arr[i, j]
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    end

    old_data_dirname = "quench_data_Jbar"
    new_data_dirname = "phase_diagram_dynamics"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram_dynamics.jld", "$(new_data_dirname)/phase_diagram_dynamics.jld")
        println("file moved")
    end

    write_log(new_data_dirname, begin_betaJ, end_betaJ, step_betaJ, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ, Λω, Λt_arr)
    export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)
end


## Tlist version
function run_phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_finecutoff(
    Tlist::Array{Float64, 1},
     begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt_arr::Array{Float64, 2}, nt::Int64,
     init_GRω::Array{ComplexF64, 3},
     mode::String = "quench", fix_beta::Bool = true, nω::Int64 = 3003)

    betaJ_num = length(Tlist)
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ
    # Sx_t_list = Array{Float64, 3}(undef, (nt + 1)÷2, betaJ_num, h_num)

    # nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    J_std = std_betaJ/(2π)
    betaJlist = 1.0 ./ (Tlist) * J_std

    if fix_beta
        println("fix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                #####################################
                # Warning: fix_beta part may not work
                #####################################
                Λt = Λt_arr[i, j]
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    else
        # println("fix beta")
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Λt = Λt_arr[i, j]
                GRωtemp = get_dynmaics_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    end

    old_data_dirname = "quench_data_Jbar"
    new_data_dirname = "phase_diagram_dynamics"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname, force=true)
        mv("phase_diagram_dynamics.jld", "$(new_data_dirname)/phase_diagram_dynamics.jld", force=true)
        println("file moved")
    end

    write_log(new_data_dirname, Tlist, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ, Λω, Λt_arr)
    export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)
end

# Using Cutoff T:  a exp(b betaJ) + bias
# Using Cutoff omega:  Λω::Float64
function run_phase_diagram_Jbar_fixβ_J_tune_init_GRω_finecutoff(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt_arr::Array{Float64, 2},
    init_GRω::Array{ComplexF64, 3}, fix_beta::Bool = true, iter_mode::Bool = false, nω::Int64 = 3003)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    ρω_centerlist = Array{ComplexF64, 4}(undef, betaJ_num, h_num, 2, 2)
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    spectral_ratio_list = Array{Float64, 2}(undef, betaJ_num, h_num)

    # nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    if fix_beta
        println("fix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                #####################################
                # Warning: fix_beta part may not work
                #####################################
                Λt = Λt_arr[i, j]
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Λω, Λt, GRω, iter_mode)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]

                spectral_ratio_list[i, j] = spectral_ratio(ρωtemp)
                println("Spectral Ratio=$(spectral_ratio_list[i, j])")

                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Λt = Λt_arr[i, j]
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, GRω, iter_mode)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]

                spectral_ratio_list[i, j] = spectral_ratio(ρωtemp)

                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    end

    old_data_dirname = "spectrum_anal_Jbar"
    new_data_dirname = "phase_diagram"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname, force=true)
        mv("phase_diagram.jld", "$(new_data_dirname)/phase_diagram.jld", force=true)
        println("file moved")
    end

    write_log(new_data_dirname, begin_betaJ, end_betaJ, step_betaJ, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ, Λω, Λt_arr)

    nt = 3003
    export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)

    # export Sxyzlist
    # Warning: Assuem we only calculate Sxyzlist that change βJ
    writedlm("$(new_data_dirname)/Sxyzlist.csv",  hcat(betaJlist, Sxyzlist[:,1,1]), ',')
    println(Sxyzlist[:,1,1])

    writedlm("$(new_data_dirname)/spectral_ratio_list.csv",  hcat(betaJlist, spectral_ratio_list[:,1]), ',')

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end

# Add spectrum ratio=( rho(0) / rho(peak) )
function run_phase_diagram_Jbar_fixβ_J_tune_init_GRω_finecutoff(
    Tlist::Array{Float64, 1},
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, Λt_arr::Array{Float64, 2},
    init_GRω::Array{ComplexF64, 3}, fix_beta::Bool = true, iter_mode::Bool = false, nω::Int64 = 3003)

    betaJ_num = length(Tlist)
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    ρω_centerlist = Array{ComplexF64, 4}(undef, betaJ_num, h_num, 2, 2)
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    spectral_ratio_list = Array{Float64, 2}(undef, betaJ_num, h_num)

    # nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω


    J_std = std_betaJ/(2π)
    betaJlist = 1.0 ./ (Tlist) * J_std

    if fix_beta
        println("fix beta is not implemented")
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        betaJlist = 1.0 ./ (Tlist) * J_std
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Λt = Λt_arr[i, j]
                ρωtemp, Sxyzlist[i, j, :], h0list[i, j], GRωtemp = get_ρω_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, GRω, iter_mode)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                center_ω = (size(ρωtemp)[1] + 1) ÷ 2
                ρω_centerlist[i, j, :, :] = ρωtemp[center_ω, :, :]

                spectral_ratio_list[i, j] = spectral_ratio(ρωtemp)
                println("Spectral Ratio=$(spectral_ratio_list[i, j])")

                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram.jld", "betaJlist",  betaJlist, "hlist", hlist, "ρω_centerlist", ρω_centerlist
            ,"Sxyzlist",Sxyzlist, "h0list", h0list)
    end

    old_data_dirname = "spectrum_anal_Jbar"
    new_data_dirname = "phase_diagram"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname, force=true)
        mv("phase_diagram.jld", "$(new_data_dirname)/phase_diagram.jld", force=true)
        println("file moved")
    end

    write_log(new_data_dirname, Tlist, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ, Λω, Λt_arr)

    nt = 3003
    export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)

    # export Sxyzlist
    # Warning: Assuem we only calculate Sxyzlist that change βJ
    writedlm("$(new_data_dirname)/Sxyzlist.csv",  hcat(betaJlist, Sxyzlist[:,1,1]), ',')
    println(Sxyzlist[:,1,1])

    writedlm("$(new_data_dirname)/spectral_ratio_list.csv",  hcat(betaJlist, spectral_ratio_list[:,1]), ',')

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end

function run_phase_diagram_dynamics_Jbar_fixβ_J_tune_init_GRω_xtoz(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
     begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
     f::Array{Float64,1}, JbarinJ::Float64, Λω::Float64, scale_Λt::Float64, nt::Int64,
     init_GRω::Array{ComplexF64, 3},
     mode::String = "quench", fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    # Sx_t_list = Array{Float64, 3}(undef, (nt + 1)÷2, betaJ_num, h_num)

    nω = size(init_GRω)[1]
    GRω_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nω, 2, 2)
    GRωtemp = Array{ComplexF64, 3}(undef, nω, 2, 2)
    GRω = init_GRω

    if fix_beta
        println("fix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                GRωtemp = get_dynmaics_Jbar_fieldxtoz(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], f, JbarinJ, beta, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    else
        # println("fix beta")
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                # Λω = scale_Λω / beta
                Λt = scale_Λt * beta
                GRωtemp = get_dynmaics_Jbar_fieldxtoz(betaJlist[i], hinJlist[j], f, JbarinJ, betaJlist[i]/J_std, Λω, Λt, nt, mode, GRω)
                if i == 1
                    GRω_smallβJ_temp[j, :,:,:] = GRωtemp
                end
                # Update initial guess
                if i == betaJ_num
                    GRω = GRω_smallβJ_temp[j,:,:,:]
                else
                    GRω = GRωtemp
                end
            end
        end
        JLD.save("phase_diagram_dynamics.jld", "betaJlist",  betaJlist, "hlist", hlist)
    end

    old_data_dirname = "quench_data_Jbar"
    new_data_dirname = "phase_diagram_dynamics"
    if isdir(new_data_dirname)
        mv(old_data_dirname, new_data_dirname)
        mv("phase_diagram_dynamics.jld", "$(new_data_dirname)/phase_diagram_dynamics.jld")
        println("file moved")
    end
end

# imaginary time equilibrium and free energy
# betaJ version doesn't add Sxiter version yet
function run_phase_diagram_Jbar_fixβ_J_imagtime_init_Gtau(
    begin_betaJ::Float64, end_betaJ::Float64, step_betaJ::Float64,
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, init_Gtau::Array{ComplexF64, 3}, fix_beta::Bool = true)

    betaJ_num = Int64(round((end_betaJ - begin_betaJ)/step_betaJ)) + 1
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    free_energylist = Array{Float64, 2}(undef, betaJ_num, h_num)
    betaJlist = collect(begin_betaJ:step_betaJ:end_betaJ)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ

    nt = size(init_Gtau)[1]
    Gt_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nt, 2, 2)
    Gttemp = Array{ComplexF64, 3}(undef, nt, 2, 2)
    Gt = init_Gtau

    if fix_beta
        println("fix beta")
        # println("\n\nfix J")
        beta = 2π
        hlist = hinJlist*std_betaJ/beta

        for j in 1:h_num
            for i in 1:betaJ_num
                #####################################
                # Warning: fix_beta part may not work
                #####################################
                Sxyzlist[i, j, :], h0list[i, j], free_energylist[i, j], Gttemp = get_Gt_imagtime_Jbar(betaJlist[i], hinJlist[j]*std_betaJ/betaJlist[i], field, f, JbarinJ, beta, Gt)
                if i == 1
                    Gt_smallβJ_temp[j, :,:,:] = Gttemp
                end

                # Update initial guess
                if i == betaJ_num
                    Gt = Gt_smallβJ_temp[j,:,:,:]
                else
                    Gt = Gttemp
                end
            end
        end
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Sxyzlist[i, j, :], h0list[i, j], free_energylist[i, j], Gttemp = get_Gt_imagtime_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Gt)
                if i == 1
                    Gt_smallβJ_temp[j, :,:,:] = Gttemp
                end

                # Update initial guess
                if i == betaJ_num
                    Gt = Gt_smallβJ_temp[j,:,:,:]
                else
                    Gt = Gttemp
                end
            end
        end
    end
    JLD.save("phase_diagram_imag.jld", "betaJlist",  betaJlist, "hlist", hlist
        ,"Sxyzlist",Sxyzlist, "h0list", h0list, "free_energylist", free_energylist)

    old_data_dirname = "Imagtime_eq_anal_Jbar"
    new_data_dirname = "phase_diagram_imag"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname, force=true)
        mv("phase_diagram_imag.jld", "$(new_data_dirname)/phase_diagram_imag.jld")
        println("file moved")
    end

    write_log(new_data_dirname, begin_betaJ, end_betaJ, step_betaJ, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ)

    # nt = 3003
    # export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)

    # export Sxyzlist
    # Warning: Assuem we only calculate Sxyzlist that change βJ
    writedlm("$(new_data_dirname)/Sxyzlist.csv",  hcat(betaJlist, Sxyzlist[:,1,1]), ',')
    println(Sxyzlist[:,1,1])

    # Warning! quick hack, assume only one column h
    writedlm("$(new_data_dirname)/free_energylist.csv",  hcat(betaJlist, free_energylist[:,1]), ',')
    println(free_energylist[:,1])

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end


# imaginary time equilibrium and free energy to have betaJlist as argument
function run_phase_diagram_Jbar_fixβ_J_imagtime_init_Gtau(
    Tlist::Array{Float64, 1},
    begin_hinJ::Float64, end_hinJ::Float64, step_hinJ::Float64, std_betaJ::Float64,
    field::String, f::Array{Float64,1}, JbarinJ::Float64, init_Gtau::Array{ComplexF64, 3},
    TorInvT::Bool, fix_beta::Bool = true, SelfEtype::String = "std", iter_mode::Bool = false, nt::Int64 = 3001)

    betaJ_num = length(Tlist)
    h_num = Int64(round((end_hinJ - begin_hinJ)/step_hinJ)) + 1
    Sxyzlist = Array{Float64, 3}(undef, betaJ_num, h_num, 4)
    h0list = Array{Float64, 2}(undef, betaJ_num, h_num)
    free_energylist = Array{Float64, 2}(undef, betaJ_num, h_num)
    hinJlist = collect(begin_hinJ:step_hinJ:end_hinJ) # before resacle to std_betaJ


    # nt = size(init_Gtau)[1]
    Gt_smallβJ_temp = Array{ComplexF64, 4}(undef, h_num, nt, 2, 2)
    Gttemp = Array{ComplexF64, 3}(undef, nt, 2, 2)
    Gt = init_Gtau

    if fix_beta
        println("fix beta mode is not implemented")
    else
        println("\n\nfix J")
        hlist = hinJlist
        J_std = std_betaJ/(2π)
        betaJlist = 1.0 ./ (Tlist) * J_std
        for j in 1:h_num
            for i in 1:betaJ_num
                beta = betaJlist[i]/J_std
                Sxyzlist[i, j, :], h0list[i, j], free_energylist[i, j], Gttemp = get_Gt_imagtime_Jbar(betaJlist[i], hinJlist[j], field, f, JbarinJ, betaJlist[i]/J_std, Gt, TorInvT, SelfEtype, iter_mode)
                if i == 1
                    Gt_smallβJ_temp[j, :,:,:] = Gttemp
                end

                # Update initial guess
                if i == betaJ_num
                    Gt = Gt_smallβJ_temp[j,:,:,:]
                else
                    Gt = Gttemp
                end
            end
        end
    end
    JLD.save("phase_diagram_imag.jld", "betaJlist",  betaJlist, "hlist", hlist
        ,"Sxyzlist",Sxyzlist, "h0list", h0list, "free_energylist", free_energylist)

    old_data_dirname = "Imagtime_eq_anal_Jbar"
    new_data_dirname = "phase_diagram_imag"
    if isdir(old_data_dirname)
        mv(old_data_dirname, new_data_dirname, force=true)
        mv("phase_diagram_imag.jld", "$(new_data_dirname)/phase_diagram_imag.jld")
        println("file moved")
    end

    write_log(new_data_dirname, Tlist, begin_hinJ, end_hinJ, step_hinJ, std_betaJ,
        field, f, JbarinJ)

    # nt = 3003
    # export_cutoff(new_data_dirname, Λω, Λt_arr, (nt-1)÷2)

    # export Sxyzlist
    # Warning: Assuem we only calculate Sxyzlist that change βJ
    writedlm("$(new_data_dirname)/Sxyzlist.csv",  hcat(betaJlist, Sxyzlist[:,1,1]), ',')
    println(Sxyzlist[:,1,1])

    # Warning! quick hack, assume only one column h
    writedlm("$(new_data_dirname)/free_energylist.csv",  hcat(betaJlist, free_energylist[:,1]), ',')
    println(free_energylist[:,1])

    # TODO, analyze the spectrum information
    # return ρωlist, Sxlist, Jbarlist
end