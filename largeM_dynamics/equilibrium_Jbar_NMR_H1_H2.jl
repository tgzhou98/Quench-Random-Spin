using LinearAlgebra
using JLD
# Warning!! Not include("miscellaneous.jl")
# This will let VSCode-Julia extension break down

function selfERealFreqSCF_Jbar(βJ::Float64, hinJ::Float64, field::String,f::Array{Float64}, JbarinJ::Float64, β::Float64,
    abscissaet::Array{Float64,1}, abscissaeω::Array{Float64,1},
    fourierTrans::Array{ComplexF64,2}, inverseTrans::Array{ComplexF64,2}, inverseTransMinust::Array{ComplexF64,2},
    nF::Array{Float64,1}, reverseNF::Array{Float64,1}, init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0)
    , Sxiter::Int64 = 0)
#   β = (2π)
    J = βJ/β
    # f = [1,1,Δ]
    b = (4/(4π*J^2))^(1/4)
    h = hinJ * J
    Jbar = JbarinJ * J
    nt = length(abscissaet)
    nω = length(abscissaeω)
    PauliMatrix = [[0 1;1 0], [0 -1im;1im 0], [1 0;0 -1]]
    HeavisideTheta(x)=x>0.0 ? 1.0 : 0.0
    GRt00 = [-1.0im*sqrt(2)*b*HeavisideTheta(t)*(1/(2*sinh(t/2))+0im)^(1/2) for t in abscissaet]
    GRω00 = fourierTrans*GRt00

    heff = h

    if init_GRω == Array{ComplexF64, 3}(undef,0,0,0)
        GRωInvfree = Array{ComplexF64, 3}(undef, nω, 2, 2)
        if field == "x"
            for i in 1:nω
                GRωInvfree[i,1,1] = abscissaeω[i] + (1e-10)im
                GRωInvfree[i,1,2] = +(heff/2)
                GRωInvfree[i,2,1] = +(heff/2)
                GRωInvfree[i,2,2] = abscissaeω[i]+(1e-10)im
            end
        elseif field == "z"
            for i in 1:nω
                GRωInvfree[i,1,1] = abscissaeω[i]+(1e-10)im+(heff/2)
                GRωInvfree[i,1,2] = 0.0
                GRωInvfree[i,2,1] = 0.0
                GRωInvfree[i,2,2] = abscissaeω[i]+(1e-10)im-(heff/2)
            end
        end

        GRω = GRωInvfree
    else
        GRω = init_GRω
    end

#   GRω= permutedims(GRω,[3,1,2]) # TODO Attention!! Different index with mathematica
  #GRω=Import[StringJoin[NotebookDirectory[],"GRω_βJ14_h0_J_x_DeltaM0_73.mx"]

  # Reshape command is anti-human (column major)
#   A = reshape(Vector(1:8), (2,2,2))
#   permutedims(A, [3, 2, 1])


  # Julia is column major, long index should put in the first
  ρω = Array{ComplexF64, 3}(undef, nω, 2, 2)
  n1 = Array{ComplexF64, 3}(undef, nω, 2, 2)
  n2 = Array{ComplexF64, 3}(undef, nω, 2, 2)
  n3 = Array{ComplexF64, 3}(undef, nω, 2, 2)
  n4 = Array{ComplexF64, 3}(undef, nω, 2, 2)
  ΣRt = Array{ComplexF64, 3}(undef, nt, 2, 2)
  ΣRω = Array{ComplexF64, 3}(undef, nω, 2, 2)
  GRωnew = Array{ComplexF64, 3}(undef, nω, 2, 2)

  dis = 100
  i = 1
  max_step = 2*10^3
#   while dis > 5*10^(-9) && i < max_step
    while dis > 10^(-6) && i < max_step
    #   while dis > 10^(-15) && i < 50

        # Choose H = - h sigma_x such that Sx > 0
        GRωInvfree = Array{ComplexF64, 3}(undef, nω, 2, 2)
        if field == "x"
            for i in 1:nω
                GRωInvfree[i,1,1] = abscissaeω[i] + (1e-10)im
                GRωInvfree[i,1,2] = +(heff/2)
                GRωInvfree[i,2,1] = +(heff/2)
                GRωInvfree[i,2,2] = abscissaeω[i]+(1e-10)im
            end
        elseif field == "z"
            for i in 1:nω
                GRωInvfree[i,1,1] = abscissaeω[i]+(1e-10)im+(heff/2)
                GRωInvfree[i,1,2] = 0.0
                GRωInvfree[i,2,1] = 0.0
                GRωInvfree[i,2,2] = abscissaeω[i]+(1e-10)im-(heff/2)
            end
        end


        for ωi in 1:nω
            ρω[ωi,:,:] = -1/(2π*1.0im)*(GRω[ωi,:,:] - adjoint(GRω[ωi,:,:]))
        end
        #Print[ρω[1:10,1,1]
        for a in 1:2
            for b in 1:2
                n1[:,a,b] = inverseTrans*(ρω[:,a,b].*reverseNF)
            end
        end

        for a in 1:2
            for b in 1:2
                n2[:,a,b] = inverseTransMinust*(ρω[:,a,b].*nF)
            end
        end

        for a in 1:2
            for b in 1:2
                n3[:,a,b] = inverseTrans*(ρω[:,a,b].*nF)
            end
        end


        for a in 1:2
            for b in 1:2
                n4[:,a,b] = inverseTransMinust*(ρω[:,a,b].*reverseNF)
            end
        end

        for ti in 1:nt
            sum_value = zeros(ComplexF64,2,2)
            for a in 1:3
                for b in 1:3
                    sum_value += (
                    (PauliMatrix[b]*n1[ti,:,:]*PauliMatrix[a])*f[a]*f[b]*tr(PauliMatrix[a]*n2[ti,:,:]*PauliMatrix[b]*n1[ti,:,:])
                    +
                    (PauliMatrix[b]*n3[ti,:,:]*PauliMatrix[a])*f[a]*f[b]*tr(PauliMatrix[a]*n4[ti,:,:]*PauliMatrix[b]*n3[ti,:,:]))
                end
            end
            ΣRt[ti,:,:] = sum_value
        end

        coeff_list = [HeavisideTheta(t) for t in abscissaet]
        for a in 1:2
            for b in 1:2
                ΣRω[:,a,b] = fourierTrans*(coeff_list.*ΣRt[:,a, b])
            end
        end

        # Convention
        # realΣRω = - J^2/4*1.0im*ΣRω[ωi,:,:]
        for ωi in 1:nω
            GRωnew[ωi,:,:] = inv((GRωInvfree[ωi,:,:] + J^2/4*1.0im*ΣRω[ωi,:,:]))
        end

        disGRω = (GRωnew - GRω)

        dis_value = 0.0
        # Warning Find another bug, 1:nω instead of nω
        for i in 1:nω
            dis_value += norm(disGRω[i,:,:]) # This p-2 norm, different from mathematica
        end
        dis = dis_value

        # Update and Mixing heff
        if Sxiter != 0

            #Prepare Thermal Equilibrium Data
            GLessert = Array{ComplexF64, 3}(undef, nt, 2, 2)
            for a in 1:2
                for b in 1:2
                    GLessert[:,a,b] = inverseTrans*(1.0im* nF.*ρω[:,a,b])
                end
            end

            # 1502, middle point
            # Sx_new = imag(GLessert[(nt+1)÷2, 1, 2])
            # Update Sx(t)
            Sxyz_new = [
            real(1/2.0*(-1.0im)*(GLessert[(nt+1)÷2, 1, 2] + GLessert[(nt+1)÷2, 2, 1])),
            real(1/2.0*(GLessert[(nt+1)÷2, 1, 2] - GLessert[(nt+1)÷2, 2, 1])),
            real(1/2.0*(-1.0im)*(GLessert[(nt+1)÷2, 1, 1] - GLessert[(nt+1)÷2, 2, 2]))
            ]
            Snorm = norm(Sxyz_new)
            append!(Sxyz_new, Snorm)
            heff_new = - Jbar * Sxyz_new[1]
            heff =  heff_new * 0.1 + heff * 0.9
        end

        GRω = 0.1*GRωnew + 0.9*GRω
        if i % 200 == 1
            println("i=$(i), dis=$(dis)")
        end
        # if i == 200
        #     println("GRω[100,1,1]:$(GRω[100,1,1]), GRω[100,1,2]:$(GRω[100,1,2]), GRω[100,2,1]:$(GRω[100,2,1]), GRω[100,2,2]:$(GRω[100,2,2])")
        #     println("ΣRω[100,1,1]:$(ΣRω[100,1,1]), ΣRω[100,1,2]:$(ΣRω[100,1,2]), ΣRω[100,2,1]:$(ΣRω[100,2,1]), ΣRω[100,2,2]:$(ΣRω[100,2,2])")
        #     println("n1[100,1,1]:$(n1[100,1,1]), n1[100,1,2]:$(n1[100,1,2]), n1[100,2,1]:$(n1[100,2,1]), n1[100,2,2]:$(n1[100,2,2])")
        # end
        i = i+1
    end
#   JLD.save("GRω_βJ$(βJ)_h$(hinJ)J_$(field)_Delta$(Δ).jld", "arr", GRω)

    # Get GRω, Sx = Im G>(t,t)

    #Prepare Thermal Equilibrium Data
    GLessert = Array{ComplexF64, 3}(undef, nt, 2, 2)
    for a in 1:2
        for b in 1:2
            GLessert[:,a,b] = inverseTrans*(1.0im* nF.*ρω[:,a,b])
        end
    end

    # 1502, middle point
    # Sx_new = imag(GLessert[(nt+1)÷2, 1, 2])
    # Update Sx(t)
    Sxyz_new = [
      real(1/2.0*(-1.0im)*(GLessert[(nt+1)÷2, 1, 2] + GLessert[(nt+1)÷2, 2, 1])),
      real(1/2.0*(GLessert[(nt+1)÷2, 1, 2] - GLessert[(nt+1)÷2, 2, 1])),
      real(1/2.0*(-1.0im)*(GLessert[(nt+1)÷2, 1, 1] - GLessert[(nt+1)÷2, 2, 2]))
    ]
    Snorm = norm(Sxyz_new)
    append!(Sxyz_new, Snorm)
    println("GLessert[(nt+1)÷2, 1, 2]", GLessert[(nt+1)÷2, 1, 2])
    # println("GLessert[(nt+1)÷2, 2, 1]", GLessert[(nt+1)÷2, 2, 1])
    println("Sxyz_new:", Sxyz_new)

    field_direct_dict = Dict("x"=>1,"y"=>2,"z"=>3)
    h0 = - (- h - f[field_direct_dict[field]] * Jbar * Sxyz_new[field_direct_dict[field]])

    if Sxiter == 0
        if i < max_step
            println("Equilibrium Converged\n")
            plot_ρω(ρω, abscissaeω, βJ, hinJ, field, f, β, 0.0)
            JLD.save("ρω.jld","arr", ρω, "GRω", GRω, "Sxyz_new", Sxyz_new, "h0", h0, "Jbar", Jbar, "omega", abscissaeω)
        else
            println("Not converged")
            # exit(1)
            plot_ρω(ρω, abscissaeω, βJ, hinJ, field, f, β, 0.0)
            JLD.save("ρω.jld","arr", ρω, "GRω", GRω, "Sxyz_new", Sxyz_new, "h0", h0, "Jbar", Jbar, "omega", abscissaeω)
        end
    else
        if i < max_step
            println("Equilibrium Converged\n")
            plot_ρω(ρω, abscissaeω, βJ, hinJ, field, f, β, 0.0, Sxiter)
            JLD.save("ρω.jld","arr", ρω, "GRω", GRω, "Sxyz_new", Sxyz_new, "h0", h0, "Jbar", Jbar, "omega", abscissaeω)
        else
            println("Not converged")
            # exit(1)
            plot_ρω(ρω, abscissaeω, βJ, hinJ, field, f, β, 0.0, Sxiter)
            JLD.save("ρω.jld","arr", ρω, "GRω", GRω, "Sxyz_new", Sxyz_new, "h0", h0, "Jbar", Jbar, "omega", abscissaeω)
        end
    end

  return ρω, Sxyz_new, h0, Jbar, GRω
end



function equilibrium_Jbar(βJ::Float64, hinJ::Float64, field::String, f::Array{Float64, 1}, JbarinJ::Float64, β::Float64, Λω::Float64, Λt::Float64, init_GRω::Array{ComplexF64, 3} = Array{ComplexF64, 3}(undef,0,0,0), iter_mode::Bool = false)

    abscissae1 = JLD.load("IntegrationData.jld", "abscissae1")
    weights1= JLD.load("IntegrationData.jld", "weights1")
    errweights1= JLD.load("IntegrationData.jld", "errweights1")

    # Λω  = 5
    nω = length(abscissae1)
    a1 = -Λω + 0.00000001
    b1 = Λω
    abscissaeω = [(b1 - a1)*abscissae1[i] + a1 for i in 1:nω]
    weightsω = [(b1 - a1)*weights1[i] for i in 1:nω]

    # Λt = 75
    nt = length(abscissae1)
    a2 = -Λt + 0.00000001
    b2 = Λt
    abscissaet = [(b2 - a2)*abscissae1[i] + a2 for i in 1:nt]
    weightst = [(b2 - a2)*weights1[i] for i in 1:nt]


    abscissaeω = vcat((abscissaeω[(1 + nω)÷2:nω] .- Λω), [0.00000000001], (Λω .+ abscissaeω[1 : (1 + nω)÷2]))
    weightsω =   vcat((weightsω[(1 + nω)÷2:nω]), [0.], (weightsω[1 : (1 + nω)÷2]))
    abscissaet = vcat((abscissaet[(1 + nt)÷2:nt] .- Λt), [0.00000000001], (Λt .+ abscissaet[1 : (1 + nt)÷2]))
    weightst =   vcat((weightst[(1 + nt)÷2 : nt]), [0.], (weightst[1 : (1 + nt)÷2]))

    weightsω = 2*Λω * weightsω/ sum(weightsω)
    weightst = 2*Λt * weightst/ sum(weightst)

    fourierTrans = Array{ComplexF64, 2}(undef, length(abscissaeω), length(abscissaet))
    for i in 1:length(abscissaeω)
        for j in 1:length(abscissaet)
            fourierTrans[i, j] = weightst[j] * exp(1.0im*abscissaet[j]*abscissaeω[i])
        end
    end


    inverseTrans= Array{ComplexF64, 2}(undef, length(abscissaet), length(abscissaeω))
    for i in 1:length(abscissaet)
        for j in 1:length(abscissaeω)
            inverseTrans[i, j] = weightsω[j] * exp(- 1.0im*abscissaet[i]*abscissaeω[j])
        end
    end


    inverseTransMinust = Array{ComplexF64, 2}(undef, length(abscissaet), length(abscissaeω))
    for i in 1:length(abscissaet)
        for j in 1:length(abscissaeω)
            inverseTransMinust[i, j] = weightsω[j] * exp(1.0im*abscissaet[i]*abscissaeω[j])
        end
    end

    nF = [1/(exp(β*ω) + 1) for ω in abscissaeω]
    reverseNF = [1/(exp(-β*ω) + 1) for ω in abscissaeω]

    Sxlist = Array{Float64, 1}(undef, 0)

    if !iter_mode
        ρω, Sxyz_new, h0, Jbar, GRω= selfERealFreqSCF_Jbar(βJ, hinJ, field, f, JbarinJ, β, abscissaet, abscissaeω, fourierTrans, inverseTrans, inverseTransMinust, nF, reverseNF, init_GRω, 0)
    else
        # Mode 1: Converge each Gt and iterative Sx in the loop
        println("Real time Sx iteration mode:")
        ρω, Sxyz_new, h0, Jbar, GRω= selfERealFreqSCF_Jbar(βJ, hinJ, field, f, JbarinJ, β, abscissaet, abscissaeω, fourierTrans, inverseTrans, inverseTransMinust, nF, reverseNF, init_GRω, 1)

        #=
        J = βJ / β
        println("Sx Iter step: 0\n")
        ρω, Sxyz_new, h0, Jbar, GRω= selfERealFreqSCF_Jbar(βJ, hinJ, field, Δ, JbarinJ, β, abscissaet, abscissaeω, fourierTrans, inverseTrans, inverseTransMinust, nF, reverseNF, init_GRω)
        Sx_old = Sxyz_new[1]
        append!(Sxlist, Sxyz_new[1])
        ρω_old = ρω
        GRω_old = GRω

        dis = 100
        i = 1
        max_step = 15
        #   while dis > 5*10^(-9) && i < max_step
        while dis > 10^(-6) && i < max_step
            hinJ_new = - JbarinJ * Sx_old
            ρω, Sxyz_new, h0, Jbar, GRω= selfERealFreqSCF_Jbar(βJ, hinJ_new, field, Δ, JbarinJ, β, abscissaet, abscissaeω, fourierTrans, inverseTrans, inverseTransMinust, nF, reverseNF, GRω_old, i)
            ρω_old = ρω
            GRω_old = GRω

            dis = abs(Sxyz_new[1] - Sx_old)
            Sx_old = 0.1* Sxyz_new[1] + 0.9 * Sx_old
            append!(Sxlist, Sx_old)
            println("Sx Iter step: $(i), dis = $(dis)\n")

            i += 1
        end
        if i < max_step
            println("Converged\n")
        else
            println("Not Converged\n")
        end
        writedlm("Sx_eq_iter.csv", Sxlist, ",")
        =#

    end

    # Plot Gt to diagnose
    nt = length(abscissaet)
    GGreatert = Array{ComplexF64, 3}(undef, nt, 2, 2)
    GLessert = Array{ComplexF64, 3}(undef, nt, 2, 2)
    for a in 1:2
        for b in 1:2
            GGreatert[:,a,b] = inverseTrans*(-1.0im*reverseNF.*ρω[:,a,b])
            GLessert[:,a,b] = inverseTrans*(1.0im* nF.*ρω[:,a,b])
        end
    end
    plot_Gt_eq(GLessert, abscissaet, βJ, hinJ, field, f, β, 0.0, 0, "equilbrium")
    JLD.save("GLessert.jld", "GLessert", GLessert)

    # Get spectrum ratio

    return ρω, abscissaeω, Sxyz_new, h0, Jbar, GRω
end



# Plot using current data
# ρω = JLD.load("ρω.jld","arr");
# plot_ρω(ρω, abscissaeω, 6.0, 0.05, "x", -0.73)