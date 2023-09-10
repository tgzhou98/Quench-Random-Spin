using MAT
using JLD
using DelimitedFiles


# Mathematica cannot import this mat file. Don't use it
function convert_all_jld2Mat()
    for (root, dirs, files) in walkdir(figdataPath)
        # println(joinpath.(root, files)) # files is a Vector{String}, can be empty
        # println(files) # files is a Vector{String}, can be empty
        absolute_path = joinpath.(root, files)
        # println(root)
        for i in 1:length(absolute_path)
            if endswith(absolute_path[i], r"imagtime_eq.jld")
                println(absolute_path[i])
                Gtau = JLD.load(absolute_path[i], "Gt")
                file = matopen(joinpath(root,"Gtau.mat"), "w")
                write(file, "Gtau", real(Gtau))
                close(file)

                tau = JLD.load(absolute_path[i], "t")
                file = matopen(joinpath(root,"tau.mat"), "w")
                write(file, "tau", tau)
                close(file)
            elseif endswith(absolute_path[i], r"ρω.jld")
                println(absolute_path[i])
                rhoomega = JLD.load(absolute_path[i], "arr")
                file = matopen(joinpath(root,"rhoomega.mat"), "w")
                write(file, "rhoomega", real(rhoomega))
                close(file)

                omega = JLD.load(absolute_path[i], "omega")
                file = matopen(joinpath(root,"omega.mat"), "w")
                write(file, "omega", omega)
                close(file)
            end
        end
    end

end


# searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

# figdataPath = "../../Rydberg-Wormhole/fig/fig4"
# key = "jld"
# convert_all_jld2Mat()
# searchdir(figdataPath, key)


# Mathematica cannot import this mat file. Don't use it
function convert_all_jld2CSV(figdataPath)
    for (root, dirs, files) in walkdir(figdataPath)
        # println(joinpath.(root, files)) # files is a Vector{String}, can be empty
        # println(files) # files is a Vector{String}, can be empty
        absolute_path = joinpath.(root, files)
        # println(root)
        for i in 1:length(absolute_path)
            if endswith(absolute_path[i], r"imagtime_eq.jld")
                println(absolute_path[i])
                Gtau = JLD.load(absolute_path[i], "Gt")
                tau = JLD.load(absolute_path[i], "t")

                writedlm(joinpath(root,"Gtau.csv"),  hcat(tau/(tau[length(tau)]), real(Gtau[:,1,1]), real(Gtau[:,1,2]), 
                                imag(Gtau[:,1,1]), imag(Gtau[:,1,2])), ',')
            elseif endswith(absolute_path[i], r"ρω.jld")
                println(absolute_path[i])
                rhoomega = JLD.load(absolute_path[i], "arr")
                omega = JLD.load(absolute_path[i], "omega")

                writedlm(joinpath(root,"rhoomega.csv"),  hcat(omega, real(rhoomega[:,1,1]), real(rhoomega[:,1,2]), 
                                real(rhoomega[:,2,1]), real(rhoomega[:,2,2])), ',')
            end
        end
    end

end


figdataPath = "../../Rydberg-Wormhole/fig/fig4"
convert_all_jld2CSV(figdataPath)

# # Pan matlab
# data_dirname = "matlab"
# if !isdir(data_dirname)
#     mkdir(data_dirname)
# end
# betaJ_string = replace("$(βJ)","."=>"_")
# target_dir = "$(data_dirname)"
# file = matopen("$(target_dir)/rhoomega.mat", "w")
# write(file, "rhoomega", ρω)
# close(file)
# file = matopen("$(target_dir)/omega.mat", "w")
# write(file, "omega", abscissaeω)
# close(file)
# file = matopen("$(target_dir)/weightsω.mat", "w")
# write(file, "weight_omega", weightsω)
# close(file)