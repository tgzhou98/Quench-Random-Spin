using LinearAlgebra
using Random
using Statistics
using ProgressMeter

PauliMatrices = [[0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]] 
Idmatrix = [1 0; 0 1] .+ 0.0im

# kron(PauliMatrices[1], PauliMatrices[1])
⊗(x, y) = kron(x, y)


############################### hamiltonian ###########################

function generate_RD_ham_Jij(N::Int64, J::Float64, Jbar::Float64, f::Array{Float64,1})

    sigmasqrt = sqrt(4*J^2/N)
    # f = [1, 1, Delta]

    H = zeros(ComplexF64, 2^N, 2^N)
    index = 1

    index = 1
    for i in 1:N
        for j in i+1:N
            # Have verified the random number is independently random at each process
            # if i == 1 && j == 2
            #     println(temp_rand)
            # end


            temp_rand = (randn() * sigmasqrt)
            #####################
            # Hack, no random
            # temp_rand = J

            for a in 1:3
                Jij_term_ops = fill(Idmatrix, N)
                Jij_term_ops[i] = PauliMatrices[a]
                Jij_term_ops[j] = PauliMatrices[a]
                # Do not introduce Jbar here. Introduce as effective h
                # temp_rand = (randn() * sigmasqrt + Jbar)

                H += (temp_rand + Jbar / N) * 1 / 4 * f[a] * foldl(⊗, Jij_term_ops)
            end

            index += 1
        end
    end

    return H
end

function generate_RD_ham_Jij_pauliXXYYZZ(N::Int64, J::Float64, Jbar::Float64)

    sigmasqrt = sqrt(4*J^2/N)
    # f = [1, 1, Delta]

    H_XXYYZZ = zeros(ComplexF64, 2^N, 2^N, 3)

    index = 1
    for i in 1:N
        for j in i+1:N
            # Have verified the random number is independently random at each process
            # if i == 1 && j == 2
            #     println(temp_rand)
            # end


            # temp_rand = (randn() * sigmasqrt)
            #####################
            # Hack, no random
            temp_rand = J

            for a in 1:3
                Jij_term_ops = fill(Idmatrix, N)
                Jij_term_ops[i] = PauliMatrices[a]
                Jij_term_ops[j] = PauliMatrices[a]
                # Do not introduce Jbar here. Introduce as effective h
                # temp_rand = (randn() * sigmasqrt + Jbar)

                H_XXYYZZ[:,:,a] += (temp_rand + Jbar / N) * 1 / 4 * foldl(⊗, Jij_term_ops)
            end

            index += 1
        end
    end

    return H_XXYYZZ
end


# Not use Jij_term_ops_arr. Directly generate H to save memory
#=
function generate_Jij_term_ops_arr(N::Int64, f::Array{Float64,1})
    Jij_term_ops_arr = zeros(ComplexF64, 2^N, 2^N, N*(N+1)÷2)
    # f = [1, 1, Delta]
    # println(f)

    index = 1
    for i in 1:N
        for j in i+1:N
            # Have verified the random number is independently random at each process
            # if i == 1 && j == 2
            #     println(temp_rand)
            # end
            for a in 1:3
                Jij_term_ops = fill(Idmatrix, N)
                Jij_term_ops[i] = PauliMatrices[a]
                Jij_term_ops[j] = PauliMatrices[a]
                # Do not introduce Jbar here. Introduce as effective h
                # temp_rand = (randn() * sigmasqrt + Jbar)
                Jij_term_ops_arr[:,:,index] += 1 / 4 * f[a] * foldl(⊗, Jij_term_ops)
            end

            index += 1
        end
    end
    
    return Jij_term_ops_arr
end

# Optimized version
function RD_ham_Jij(N::Int64, J::Float64, Jbar::Float64, f::Array{Float64,1}, Jij_term_ops_arr::Array{ComplexF64, 3})

    sigmasqrt = sqrt(4*J^2/N)
    # f = [1, 1, Delta]

    H = zeros(ComplexF64, 2^N, 2^N)
    index = 1
    for i in 1:N
        for j in i+1:N
            # temp_rand = (randn() * sigmasqrt)
            #####################
            # Hack, no random
            temp_rand = J

            # Have verified the random number is independently random at each process
            # if i == 1 && j == 2
            #     println(temp_rand)
            # end
            H += (temp_rand + Jbar / N) * Jij_term_ops_arr[:,:,index]
            index += 1
        end
    end

    return H
end
=#

function RD_ham_h(N::Int64, J::Float64, hinJ::Float64)

    h = hinJ * J
    H = zeros(ComplexF64, 2^N, 2^N)

    a = 1
    for i in 1:N
        hi_term_ops = fill(Idmatrix, N)
        hi_term_ops[i] = PauliMatrices[a]
        H -= h / 2 * foldl(⊗, hi_term_ops)
    end

    return H
end

############################# finite temperature calculation ############################

###################################### Equilibrium ######################################
function equilibrium_Sx(disorder_num::Int64, beta::Float64, N::Int64, J::Float64, 
        f::Array{Float64,1}, hinJ::Float64, Jbar::Float64, Jij_term_ops_arr::Array{ComplexF64, 3})
    Sx_list = zeros(disorder_num)
    # @showprogress 1 "Computing disorder" for dis_index in 1:disorder_num

    a = 1
    Sx_op = zeros(ComplexF64, 2^N, 2^N)
    for i in 1:N
        hi_term_ops = fill(Idmatrix, N)
        hi_term_ops[i] = PauliMatrices[a]
        Sx_op += 1 / 2 * foldl(⊗, hi_term_ops)
    end
    # Hamh0 = RD_ham_h(N, J, hinJ)
    Hamh0 = (-hinJ * J) * Sx_op
    for dis_index in 1:disorder_num
        # if dis_index % 20 == 1
        #     println("$(dis_index)/$(disorder_num)")
        # end
        HamJij = RD_ham_Jij(N, J, Jbar, f, Jij_term_ops_arr)
        Ham1 = HamJij + Hamh0
        rho = exp(- beta*Ham1)
        rho = rho / tr(rho)

        Sx_list[dis_index] = real(tr(rho*Sx_op)) / N
    end
    return mean(Sx_list)
end



