using JLD

function Cutoffomega()
    data_dirname = "../data/analyze_dynamics/phase_diagram_dynamics_Jbar_alpha_M0_73_high_accu_fixJ_init_GR/Jbar0_0"
    hinJ = 0.2
    cutoff_omega = []
    for βJ in collect(4.0:14.0)
        betaJ_string = replace("$(round(βJ, digits=2))","."=>"_")
        hinJ_string = replace("$(round(hinJ, digits=3))","."=>"_")
        target_dir = "$(data_dirname)/betaJ$(betaJ_string)_hinJ$(hinJ_string)"
        
        rhoomega_data_path = "$(target_dir)/ρω.jld"
        rhoomega_data = JLD.load(rhoomega_data_path, "arr")
        abscissaeω = JLD.load(rhoomega_data_path, "omega")
        # Gt_data = JLD.load(rhoomega_data_path, "Gt")

        nω = size(rhoomega_data)[1]

        final_index = 1
        for i in 1:nω
            if abs(rhoomega_data[i,1,1]) > 10^(-2)
                final_index = i
            end
        end

        append!(cutoff_omega, abscissaeω[final_index])
        println(final_index)
        println(abscissaeω[nω])
    end

    return cutoff_omega
end

println(Cutoffomega())