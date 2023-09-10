using LinearAlgebra
using JLD
using DataFrames

include("SY_dynamics.jl")
PauliMatrix = [[0 1;1 0], [0 -1im;1im 0], [1 0;0 -1]]

# BUG FIXED
# Self Energy must add a extra minus sign

# auxilary function to update energy in spin sum


function self_energy_update(i::Int64, j::Int64, J::Float64, f::Array{Float64, 1},
  temp_selfE_Greater::Array{ComplexF64, 2}, temp_selfE_Lesser::Array{ComplexF64, 2},
  GGreaterDynamicsttp::Array{ComplexF64, 4}, GLesserDynamicsttp::Array{ComplexF64, 4},
  ΣGreaterDynamicsttp::Array{ComplexF64, 4}, ΣLesserDynamicsttp::Array{ComplexF64, 4},
  )
    fill!(temp_selfE_Greater, zero(ComplexF64))
    fill!(temp_selfE_Lesser, zero(ComplexF64))
    for a in 1:3
      for b in 1:3
        temp_selfE_Greater += (J^2/4)*((PauliMatrix[b]*GGreaterDynamicsttp[i, j, :, :]*PauliMatrix[a])*f[a]*f[b]
          *tr(PauliMatrix[a]*GLesserDynamicsttp[j, i, :, :]*PauliMatrix[b]*GGreaterDynamicsttp[i, j, :, :])
          )
        temp_selfE_Lesser += (J^2/4)*((PauliMatrix[b]*GLesserDynamicsttp[i, j, :, :]*PauliMatrix[a])*f[a]*f[b]
        *tr(PauliMatrix[a]*GGreaterDynamicsttp[j, i, :, :]*PauliMatrix[b]*GLesserDynamicsttp[i, j, :, :])
        )
      end
    end
    ΣGreaterDynamicsttp[i,j,:,:] = temp_selfE_Greater
    ΣLesserDynamicsttp[i,j,:,:]  = temp_selfE_Lesser
end

function self_energy_update(i::Int64, j::Int64, J::Float64, f::Array{Float64, 1}, type::String,
  temp_selfE_Greater::Array{ComplexF64, 2}, temp_selfE_Lesser::Array{ComplexF64, 2},
  GGreaterDynamicsttp::Array{ComplexF64, 4}, GLesserDynamicsttp::Array{ComplexF64, 4},
  ΣGreaterDynamicsttp::Array{ComplexF64, 4}, ΣLesserDynamicsttp::Array{ComplexF64, 4},
  ΣRetardedDynamicsttp::Array{ComplexF64, 4}, GRetardedDynamicsttp::Array{ComplexF64, 4},
  ΣAdavancedDynamicsttp::Array{ComplexF64, 4}, GAdavancedDynamicsttp::Array{ComplexF64, 4},
  )
    fill!(temp_selfE_Greater, zero(ComplexF64))
    fill!(temp_selfE_Lesser, zero(ComplexF64))
    for a in 1:3
      for b in 1:3
        temp_selfE_Greater += (J^2/4)*((PauliMatrix[b]*GGreaterDynamicsttp[i, j, :, :]*PauliMatrix[a])*f[a]*f[b]
          *tr(PauliMatrix[a]*GLesserDynamicsttp[j, i, :, :]*PauliMatrix[b]*GGreaterDynamicsttp[i, j, :, :])
          )
        temp_selfE_Lesser += (J^2/4)*((PauliMatrix[b]*GLesserDynamicsttp[i, j, :, :]*PauliMatrix[a])*f[a]*f[b]
        *tr(PauliMatrix[a]*GGreaterDynamicsttp[j, i, :, :]*PauliMatrix[b]*GLesserDynamicsttp[i, j, :, :])
        )
      end
    end
    if i == 1001 && j == 1000
      println(temp_selfE_Greater)
    end

    ΣGreaterDynamicsttp[i,j,:,:] = temp_selfE_Greater
    ΣLesserDynamicsttp[i,j,:,:]  = temp_selfE_Lesser

    if type == "down"
      ΣRetardedDynamicsttp[i, j,:,:] =  (ΣGreaterDynamicsttp[i,j,:,:] - ΣLesserDynamicsttp[i,j,:,:])
      GRetardedDynamicsttp[i, j,:,:] = (GGreaterDynamicsttp[i,j,:,:] - GLesserDynamicsttp[i,j,:,:])
    elseif type == "right"
      ΣAdavancedDynamicsttp[i, j,:,:] =  (ΣLesserDynamicsttp[i,j,:,:] - ΣGreaterDynamicsttp[i,j,:,:])
      GAdavancedDynamicsttp[i, j,:,:] = (GLesserDynamicsttp[i,j,:,:] - GGreaterDynamicsttp[i,j,:,:] )
    elseif type == "corner"
      ΣRetardedDynamicsttp[i, j,:,:] =  (ΣGreaterDynamicsttp[i,j,:,:] - ΣLesserDynamicsttp[i,j,:,:])*0.5
      GRetardedDynamicsttp[i, j,:,:] = (GGreaterDynamicsttp[i,j,:,:] - GLesserDynamicsttp[i,j,:,:])*0.5

      ΣAdavancedDynamicsttp[i, j,:,:] =  (ΣLesserDynamicsttp[i,j,:,:] - ΣGreaterDynamicsttp[i,j,:,:] )*0.5
      GAdavancedDynamicsttp[i, j,:,:] = (GLesserDynamicsttp[i,j,:,:] - GGreaterDynamicsttp[i,j,:,:] )*0.5
    end
end

function baym_kadanoff_Jbar(J::Float64, h::Float64, field::String, Sxyz_old::Array{Float64, 1}
      ,Jbar::Float64, nt::Int64, Δt::Float64, f::Array{Float64,1}
      ,GGreaterDynamicsttp::Array{ComplexF64,4}
      ,GLesserDynamicsttp::Array{ComplexF64,4}
      ,ΣGreaterDynamicsttp::Array{ComplexF64,4}
      ,ΣLesserDynamicsttp::Array{ComplexF64,4}
      ,GRetardedDynamicsttp::Array{ComplexF64,4}
      ,GAdavancedDynamicsttp::Array{ComplexF64,4}
      ,ΣRetardedDynamicsttp::Array{ComplexF64,4}
      ,ΣAdavancedDynamicsttp::Array{ComplexF64,4}
      ,mode::String
      ,verbose::Bool
)
    #baym_kadanoff
    # order 2 Runge-Kutta iteration doesn't work. It couple Greater and Lesser Green function.
    # Predictor Correlator won't work 

  criterion = abs(1e-5 * norm(Sxyz_old))

  # Choose H0 such that Sx > 0
  H0 = zeros(ComplexF64, 2, 2)
  # Update H0 in the beginning of the iteration
  for i in 1:3
    H0 += f[i]*PauliMatrix[i]*Sxyz_old[i]*Jbar / 2.0
  end
  if field == "x"
    H0 += -PauliMatrix[1] * h / 2
  elseif field == "z"
    H0 += -PauliMatrix[3] * h / 2
  end

  cutoff_length = (nt - 1)÷2
  # Deriv_G1_Greater iter*save_deriv*2*2
  # Deriv_G_diag_Greater save_deriv*2*2
  Deriv_G1_Greater = Array{ComplexF64, 4}(undef, cutoff_length,2, 2, 2)
  Deriv_G2_Greater = Array{ComplexF64, 4}(undef, cutoff_length,2, 2, 2)
  Deriv_G1_Lesser = Array{ComplexF64, 4}(undef, cutoff_length,2, 2, 2)
  Deriv_G2_Lesser = Array{ComplexF64, 4}(undef, cutoff_length,2, 2, 2)
  Deriv_G_diag_Greater = Array{ComplexF64, 3}(undef,2, 2, 2)
  Deriv_G_diag_Lesser  = Array{ComplexF64, 3}(undef,2, 2, 2)

  ######################################################################################################
  # Initial Deriv
  ######################################################################################################
  k = (nt + 1)÷2
  println(k)
  @time begin
  GR_pas = GRetardedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
  GA_pas = GAdavancedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
  GGreater_pas = GGreaterDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
  GLesser_pas = GLesserDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]

  ΣR_pas = ΣRetardedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
  ΣA_pas =ΣAdavancedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
  ΣGreater_pas =ΣGreaterDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:] 
  ΣLesser_pas =ΣLesserDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:] 

  relative_index = 1

  # Tranpose is a lazy wrapper type of transpose

  # Initial temp_sum
  temp_sum = zeros(ComplexF64, cutoff_length)

  # G1 >
  G1_Greater_deriv(temp_sum, ΣR_pas, GGreater_pas, ΣGreater_pas, GA_pas,
    Deriv_G1_Greater, H0, Δt, relative_index, cutoff_length)

  # G1 <
  G1_Lesser_deriv(temp_sum, ΣR_pas, GLesser_pas, ΣLesser_pas, GA_pas,
    Deriv_G1_Lesser, H0, Δt, relative_index, cutoff_length)

  # G2 >
  G2_Greater_deriv(temp_sum, GR_pas, ΣGreater_pas, GGreater_pas, ΣA_pas,
    Deriv_G2_Greater, H0, Δt, relative_index, cutoff_length)

  # G2 <
  G2_Lesser_deriv(temp_sum, GR_pas, ΣLesser_pas, GLesser_pas, ΣA_pas,
    Deriv_G2_Lesser, H0, Δt, relative_index, cutoff_length)

  # Diagnoal G(t,t) evolution
  # TODO: !!!! \partial t_1 + \partial t_2 = \frac{d}{dt}
  # i \frac{d}{dt} G(t,t) = I_1(t,t) - I_2(t,t)
  temp_sum = zero(ComplexF64)
  # Diagnoal
  G_diagonal_Greater_deriv(temp_sum, ΣR_pas, GGreater_pas, ΣGreater_pas, GA_pas, GR_pas, ΣA_pas,
      Deriv_G_diag_Greater, H0, Δt, relative_index, cutoff_length)

  G_diagonal_Lesser_deriv(temp_sum, ΣR_pas, GLesser_pas, ΣLesser_pas, GA_pas, GR_pas, ΣA_pas,
      Deriv_G_diag_Lesser, H0, Δt, relative_index, cutoff_length)

  end

  k_save = 0
  
  ######################################################################################################
  # Start Iteration
  ######################################################################################################
  for k in (nt + 1)÷2:nt - 1
    H0 = zeros(ComplexF64, 2, 2)
    # Update H0 in the beginning of the iteration
    for i in 1:3
      H0 += f[i]*PauliMatrix[i]*Sxyz_old[i]*Jbar / 2.0
    end
    if field == "x"
      H0 += -PauliMatrix[1] * h / 2
    elseif field == "z"
      H0 += -PauliMatrix[3] * h / 2
    end

    GR_pas = GRetardedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
    GA_pas = GAdavancedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
    GGreater_pas = GGreaterDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
    GLesser_pas = GLesserDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]

    ΣR_pas = ΣRetardedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
    ΣA_pas =ΣAdavancedDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:]
    ΣGreater_pas =ΣGreaterDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:] 
    ΣLesser_pas =ΣLesserDynamicsttp[k-cutoff_length+1:k,k-cutoff_length+1:k,:,:] 

    ######################################################################################################
    # Update the new deriviative
    ######################################################################################################

    ############################################################################
    ## relative_index = 1 Update the end
    ############################################################################
    relative_index = 1
      # Initial temp_sum
      temp_sum = zeros(ComplexF64, cutoff_length)

      # G1 >
      G1_Greater_deriv(temp_sum, ΣR_pas, GGreater_pas, ΣGreater_pas, GA_pas,
        Deriv_G1_Greater, H0, Δt, relative_index, cutoff_length,true)

      # G1 <
      G1_Lesser_deriv(temp_sum, ΣR_pas, GLesser_pas, ΣLesser_pas, GA_pas,
        Deriv_G1_Lesser, H0, Δt, relative_index, cutoff_length,true)

      # G2 >
      G2_Greater_deriv(temp_sum, GR_pas, ΣGreater_pas, GGreater_pas, ΣA_pas,
        Deriv_G2_Greater, H0, Δt, relative_index, cutoff_length,true)

      # G2 <
      G2_Lesser_deriv(temp_sum, GR_pas, ΣLesser_pas, GLesser_pas, ΣA_pas,
        Deriv_G2_Lesser, H0, Δt, relative_index, cutoff_length,true)


    ############################################################################
    ## relative_index = 2
    ############################################################################
    relative_index = 2

      # Tranpose is a lazy wrapper type of transpose
      # Initial temp_sum
      temp_sum = zeros(ComplexF64, cutoff_length)

      # G1 >
      G1_Greater_deriv(temp_sum, ΣR_pas, GGreater_pas, ΣGreater_pas, GA_pas,
        Deriv_G1_Greater, H0, Δt, relative_index, cutoff_length)

      # G1 <
      G1_Lesser_deriv(temp_sum, ΣR_pas, GLesser_pas, ΣLesser_pas, GA_pas,
        Deriv_G1_Lesser, H0, Δt, relative_index, cutoff_length)

      # G2 >
      G2_Greater_deriv(temp_sum, GR_pas, ΣGreater_pas, GGreater_pas, ΣA_pas,
        Deriv_G2_Greater, H0, Δt, relative_index, cutoff_length)

      # G2 <
      G2_Lesser_deriv(temp_sum, GR_pas, ΣLesser_pas, GLesser_pas, ΣA_pas,
        Deriv_G2_Lesser, H0, Δt, relative_index, cutoff_length)

      # Diagnoal G(t,t) evolution
      # TODO: !!!! \partial t_1 + \partial t_2 = \frac{d}{dt}
      # i \frac{d}{dt} G(t,t) = I_1(t,t) - I_2(t,t)
      temp_sum = zero(ComplexF64)
      # Diagnoal
      G_diagonal_Greater_deriv(temp_sum, ΣR_pas, GGreater_pas, ΣGreater_pas, GA_pas, GR_pas, ΣA_pas,
          Deriv_G_diag_Greater, H0, Δt, relative_index, cutoff_length)

      G_diagonal_Lesser_deriv(temp_sum, ΣR_pas, GLesser_pas, ΣLesser_pas, GA_pas, GR_pas, ΣA_pas,
          Deriv_G_diag_Lesser, H0, Δt, relative_index, cutoff_length)


    if verbose
      println("upup dG1(k, k-1)   2   save before            ", (Deriv_G1_Greater[(nt-1)÷2-1,2,1,1])*Δt)
      println("upup dG1(k, k-1)   1   save before            ", (Deriv_G1_Greater[(nt-1)÷2-1,1,1,1])*Δt)
    end

    ######################################################################################################
    #Start Down line
    ######################################################################################################
    #Update Greater DownList
    GGreaterDynamicsttp[k + 1, k - (nt - 1)÷2 + 1 : k,:,:] = (3/2 * Deriv_G1_Greater[:,2,:,:] - 1/2 * Deriv_G1_Greater[:,1,:,:]) * Δt + GGreaterDynamicsttp[k, k - (nt - 1)÷2 + 1 : k,:,:]
    #Update Lesser DownList
    GLesserDynamicsttp[k + 1, k - (nt - 1)÷2 + 1 : k,:,:] = (3/2 * Deriv_G1_Lesser[:,2,:,:] - 1/2 * Deriv_G1_Lesser[:,1,:,:]) * Δt + GLesserDynamicsttp[k, k - (nt - 1)÷2 + 1 : k,:,:]


    ######################################################################################################
    #Start Right line
    ######################################################################################################
    #Update Greater DownList
    GGreaterDynamicsttp[k - (nt - 1)÷2 + 1 : k, k + 1,:,:] = (3/2 * Deriv_G2_Greater[:,2,:,:] - 1/2 * Deriv_G2_Greater[:,1,:,:]) * Δt + GGreaterDynamicsttp[k - (nt - 1)÷2 + 1 : k, k,:,:]
    #Update Lesser DownList
    GLesserDynamicsttp[k - (nt - 1)÷2 + 1 : k, k + 1,:,:] = (3/2 * Deriv_G2_Lesser[:,2,:,:] - 1/2 * Deriv_G2_Lesser[:,1,:,:]) * Δt + GLesserDynamicsttp[k - (nt - 1)÷2 + 1 : k, k,:,:]


    ######################################################################################################
    # Diagnoal Term
    ######################################################################################################
    GGreaterDynamicsttp[k + 1, k + 1,:,:] = GGreaterDynamicsttp[k, k,:,:] + (3/2*Deriv_G_diag_Greater[2,:,:] - 1/2*Deriv_G_diag_Greater[1,:,:]) * Δt
    GLesserDynamicsttp[k + 1, k + 1,:,:] =  GLesserDynamicsttp[k, k,:,:] + (3/2*Deriv_G_diag_Lesser[2,:,:] - 1/2*Deriv_G_diag_Lesser[1,:,:]) * Δt

    # Hard code diagonal term 

    # GGreaterDynamicsttp[k + 1, k + 1,:,:] = GGreaterDynamicsttp[1, 1,:,:] 
    # GLesserDynamicsttp[k + 1, k + 1,:,:] = GLesserDynamicsttp[1, 1,:,:] 

    if verbose

    # println("upup dG1(k, k-2:k)                   ", Deriv_G1_Greater[(nt-1)÷2-2:(nt-1)÷2,1,1]*Δt)
    # println("upup G>pas(k,i) - G>pas(k,i-1)       ", (GGreaterDynamicsttp[k, k,1,1] - GGreaterDynamicsttp[k, k-1,1,1]))
    # println("upup dG1(k, k-2:k)                   ", Deriv_G1_Greater[(nt-1)÷2-2:(nt-1)÷2,1,1]*Δt)

    println("upup dG1(k, k-1)                     ", (3/2*Deriv_G1_Greater[(nt-1)÷2-1,2,1,1] - 1/2*Deriv_G1_Greater[(nt-1)÷2-1,1,1,1])*Δt)
    println("upup G>pas(k+1,k)                    ",    GGreaterDynamicsttp[k + 1, k,1,1])
    println("upup G>pas(k+1,k) - G>pas(k+1,k-1)   ",(GGreaterDynamicsttp[k + 1, k,1,1] - GGreaterDynamicsttp[k + 1, k-1,1,1]))

    # println("upup dG1(k, k-2)                     ", Deriv_G1_Greater[(nt-1)÷2-2,1,1]*Δt)
    # println("upup G>pas(k+1,k-1)                  ",    GGreaterDynamicsttp[k + 1, k-1,1,1])
    # println("upup G>pas(k+1,k-1) - G>pas(k+1,k-2) ",(GGreaterDynamicsttp[k + 1, k-1,1,1] - GGreaterDynamicsttp[k + 1, k-2,1,1]))

    println("upup dG2(k-1, k)                     ", (3/2*Deriv_G2_Greater[(nt-1)÷2-1,2,1,1] - 1/2*Deriv_G2_Greater[(nt-1)÷2-1,1,1,1])*Δt)
    println("upup G>pas(k, k+1)                   ",    GGreaterDynamicsttp[k, k + 1,1,1])
    println("upup G>pas(k, k+1) - G>pas(k-1, k+1) ",(GGreaterDynamicsttp[k, k + 1,1,1] - GGreaterDynamicsttp[k-1, k + 1,1,1]))

    # println("updown dG1(k, k-2:k)                   ", Deriv_G1_Greater[(nt-1)÷2-2:(nt-1)÷2,1,2]*Δt)
    # println("updown G>pas(k,i) - G>pas(k,i-1)       ", (GGreaterDynamicsttp[k, k,1,2] - GGreaterDynamicsttp[k, k-1,1,2]))

    println("updown dG1(k, k-1)                     ", (3/2*Deriv_G1_Greater[(nt-1)÷2-1,2,1,2] - 1/2*Deriv_G1_Greater[(nt-1)÷2-1,1,1,2])*Δt)
    println("updown G>pas(k+1,k)                    ",    GGreaterDynamicsttp[k + 1, k,1,2])
    println("updown G>pas(k+1,k) - G>pas(k+1,k-1)   ",(GGreaterDynamicsttp[k + 1, k,1,2] - GGreaterDynamicsttp[k + 1, k-1,1,2]))

    # println("updown dG1(k, k-2)                     ", Deriv_G1_Greater[(nt-1)÷2-2,1,2]*Δt)
    # println("updown G>pas(k+1,k-1)                  ",    GGreaterDynamicsttp[k + 1, k-1,1,2])
    # println("updown G>pas(k+1,k-1) - G>pas(k+1,k-2) ",(GGreaterDynamicsttp[k + 1, k-1,1,2] - GGreaterDynamicsttp[k + 1, k-2,1,2]))

    println("updown dG2(k-1, k)                     ", (3/2.0*Deriv_G2_Greater[(nt-1)÷2-1,2,1,2] - 1/2.0*Deriv_G2_Greater[(nt-1)÷2-1,1,1,2])*Δt)
    println("updown G>pas(k, k+1)                   ",    GGreaterDynamicsttp[k, k + 1,1,2])
    println("updown G>pas(k, k+1) - G>pas(k-1, k+1) ",(GGreaterDynamicsttp[k, k + 1,1,2] - GGreaterDynamicsttp[k-1, k + 1,1,2]))

    end

    # println("upup dG2(k-1, k)                     ", Deriv_G2_Greater[(nt-1)÷2,1,1]*Δt)
    # println("upup G>pas(k, k+1)                 ",    GGreaterDynamicsttp[k, k + 1,1,1])
    # println("upup G>pas(k,k+1) - G>pas(k-1,k+1)    ",(GGreaterDynamicsttp[k, k + 1,1,1] - GGreaterDynamicsttp[k-1, k + 1,1,1]))

    if verbose
      println("Deriv_G_diag_Greater", (3/2*Deriv_G_diag_Greater[2,1,1] - 1/2*Deriv_G_diag_Greater[1,1,1]))
      println("Deriv_G_diag_Lesser", (3/2*Deriv_G_diag_Lesser[2,1,1] - 1/2*Deriv_G_diag_Lesser[1,1,1]))
    end


    # if verbose
    #   println("upup G>(k+1,k+1)                 ",    GGreaterDynamicsttp[k + 1, k+1,1,1])
    #   println("upup G<(k+1,k+1)                 ",    GLesserDynamicsttp[k + 1, k+1,1,1])
    # end
    

    # Initial temp memory for spin sum
    temp_selfE_Greater = zeros(ComplexF64, 2, 2)
    temp_selfE_Lesser = zeros(ComplexF64, 2, 2)
    
    # Update Self Energy(>, <, R, A) and Green function(R,A)
    i = k+1
    for j in k - (nt - 1)÷2 + 1 : k
      self_energy_update(i, j, J, f, "down", temp_selfE_Greater, temp_selfE_Lesser,
        GGreaterDynamicsttp, GLesserDynamicsttp, ΣGreaterDynamicsttp, ΣLesserDynamicsttp,
        ΣRetardedDynamicsttp, GRetardedDynamicsttp, ΣAdavancedDynamicsttp, GAdavancedDynamicsttp)
    end

    j = k+1
    for i in k - (nt - 1)÷2 + 1 : k
      self_energy_update(i, j, J, f,"right", temp_selfE_Greater, temp_selfE_Lesser,
        GGreaterDynamicsttp, GLesserDynamicsttp, ΣGreaterDynamicsttp, ΣLesserDynamicsttp,
        ΣRetardedDynamicsttp, GRetardedDynamicsttp, ΣAdavancedDynamicsttp, GAdavancedDynamicsttp)
    end

    i = k+1
    j = k+1
    self_energy_update(i, j, J, f,"corner", temp_selfE_Greater, temp_selfE_Lesser,
      GGreaterDynamicsttp, GLesserDynamicsttp, ΣGreaterDynamicsttp, ΣLesserDynamicsttp,
      ΣRetardedDynamicsttp, GRetardedDynamicsttp, ΣAdavancedDynamicsttp, GAdavancedDynamicsttp)

    ######################################################################################################
    #TODO Update Derivative index
    # Attention!! Not simply assign upper line to down line. Only cutoff_length - 1?? can be reuse
    ######################################################################################################
    Deriv_G1_Greater[1:cutoff_length-1,1,:,:]     =  Deriv_G1_Greater[2:cutoff_length,2,:,:]    
    Deriv_G2_Greater[1:cutoff_length-1,1,:,:]     =  Deriv_G2_Greater[2:cutoff_length,2,:,:]    
    Deriv_G1_Lesser[1:cutoff_length-1,1,:,:]      =  Deriv_G1_Lesser[2:cutoff_length,2,:,:]     
    Deriv_G2_Lesser[1:cutoff_length-1,1,:,:]      =  Deriv_G2_Lesser[2:cutoff_length,2,:,:]     
    Deriv_G_diag_Greater[1,:,:] =  Deriv_G_diag_Greater[2,:,:]
    Deriv_G_diag_Lesser[1,:,:]  =  Deriv_G_diag_Lesser[2,:,:] 


    if verbose
      println("upup dG1(k, k-1)   2   save after            ", (Deriv_G1_Greater[(nt-1)÷2-1,2,1,1])*Δt)
      println("upup dG1(k, k-1)   1   save after            ", (Deriv_G1_Greater[(nt-1)÷2-1,1,1,1])*Δt)
    end

    # Update Sx(t)
    Sxyz_new = [
      real(1/2.0*(-1.0im)*(GLesserDynamicsttp[k+1, k+1, 1, 2] + GLesserDynamicsttp[k+1, k+1, 2, 1])),
      real(1/2.0*(GLesserDynamicsttp[k+1, k+1, 1, 2] - GLesserDynamicsttp[k+1, k+1, 2, 1])),
      real(1/2.0*(-1.0im)*(GLesserDynamicsttp[k+1, k+1, 1, 1] - GLesserDynamicsttp[k+1, k+1, 2, 2]))
    ]
    Snorm = norm(Sxyz_new)
    append!(Sxyz_new, Snorm)
    Sxyz_old = Sxyz_new

    if k % 200 == 0
      println("Step:$k, Sx:$(Sxyz_old[1]), Sy:$(Sxyz_old[2]), Sz:$(Sxyz_old[3]), |S|:$(Sxyz_old[4])")
    end 
    k_save = k

    if mode == "quench"
      if (Sxyz_old[4] < criterion && (k - (nt + 1)÷2) / cutoff_length > 2 / 5) || (k - (nt + 1)÷2) / cutoff_length > 3 / 5 || Sxyz_old[4] > 2.0 
      # if ((k - (nt + 1)÷2) >= cutoff_length - 4) || Sxyz_old[4] > 2.0
      # if (Sxyz_old[4] > 2.0 && (Sxyz_old[4]) < criterion && (k - (nt + 1)÷2) / cutoff_length > 1 / 5)
        # JLD.save("Dynamics.jld", "GGreaterDynamicsttp",  GGreaterDynamicsttp, "GLesserDynamicsttp", GLesserDynamicsttp)
        break
      end
    elseif mode == "TI"
      if (k - (nt + 1)÷2) / cutoff_length > 1 / 5
        # JLD.save("Dynamics.jld", "GGreaterDynamicsttp",  GGreaterDynamicsttp, "GLesserDynamicsttp", GLesserDynamicsttp)
        break
      end
    end
  end
  return k_save
    
end

function dynamics_init_Jbar(βJ::Float64, field::String, f::Array{Float64,1}, β::Float64, ρω::Array{ComplexF64, 3}, Λω::Float64, Λt::Float64, nt::Int64 = 2001)

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
  a2=-Λt+0.00000001;
  b2= Λt;
  # nt = 2001;
  equaldTabscissae = collect(0:1/(nt-1):1)
  Δt=Λt/((nt-1)/2);
  equaldTabscissaet=[(b2-a2)*equaldTabscissae[i]+a2 for i in 1:nt];
  # equaldTweightst=Λt*equaldTweights/sum(equaldTweights);

  abscissaeω = vcat((abscissaeω[(1 + nω)÷2:nω] .- Λω), [0.00000000001], (Λω .+ abscissaeω[1 : (1 + nω)÷2]))
  weightsω =   vcat((weightsω[(1 + nω)÷2:nω]), [0.], (weightsω[1 : (1 + nω)÷2]))

  weightsω = 2*Λω * (weightsω/ sum(weightsω))

  invEqualdTTrans= Array{ComplexF64, 2}(undef, length(equaldTabscissaet), length(abscissaeω))
  for i in 1:length(equaldTabscissaet)
      for j in 1:length(abscissaeω)
        invEqualdTTrans[i, j] = weightsω[j] * exp(- 1.0im*equaldTabscissaet[i]*abscissaeω[j])
      end
  end

  nF = [1/(exp(β*ω) + 1) for ω in abscissaeω]
  reverseNF = [1/(exp(-β*ω) + 1) for ω in abscissaeω]

  # β = (2π)
  J = βJ/β
  # h = hinJ * J
  # f = [1,1,Δ]
  b = (4/(4π*J^2))^(1/4)
  # if field == "x"
  #   H0 = [0 h/2; h/2 0] .+ 0im
  # elseif field == "z"
  #   H0 = [h/2 0; 0 -h/2] .+ 0im
  # end

  #Prepare Thermal Equilibrium Data
  GGreatert = Array{ComplexF64, 3}(undef, nt, 2, 2)
  GLessert = Array{ComplexF64, 3}(undef, nt, 2, 2)
  for a in 1:2
      for b in 1:2
          GGreatert[:,a,b] = invEqualdTTrans*(-1.0im*reverseNF.*ρω[:,a,b])
          GLessert[:,a,b] = invEqualdTTrans*(1.0im* nF.*ρω[:,a,b])
      end
  end
  showstep = 1
  tempnomega = length(abscissaeω)
  println("###################################################")
  println("Start printing initial ρω solution")
  display(ρω[(tempnomega+1)÷2-showstep:(tempnomega+1)÷2+showstep,:,:])
  display(nF[(tempnomega+1)÷2-showstep:(tempnomega+1)÷2+showstep,:,:])
  display(reverseNF[(tempnomega+1)÷2-showstep:(tempnomega+1)÷2+showstep,:,:])
  println("###################################################")

  println("GLessert in dynamics")
  display(GLessert[(nt+1)÷2-showstep:(nt+1)÷2+showstep, 1, 2])
  display(GGreatert[(nt+1)÷2-showstep:(nt+1)÷2+showstep, 1, 2])

  #Start Filling initial condition
  GGreaterDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 
  GLesserDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 
  for i in 1:(nt+1)÷2
    for j in 1:(nt+1)÷2
      GGreaterDynamicsttp[i, j, :, :] = GGreatert[i + (nt + 1)÷2 - j, :, :]
      GLesserDynamicsttp[i, j, :, :] = GLessert[i + (nt + 1)÷2 - j, :, :]
    end
  end
  println("G< up down", GLessert[(nt + 1)÷2, 1, 2])
  
  #Start initial self energy
  ΣGreaterDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 
  ΣLesserDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 


  # Initial temp memory for spin sum
  temp_selfE_Greater = zeros(ComplexF64, 2, 2)
  temp_selfE_Lesser = zeros(ComplexF64, 2, 2)

  for i in 1:(nt+1)÷2
    for j in 1:(nt+1)÷2
      self_energy_update(i, j, J, f, temp_selfE_Greater, temp_selfE_Lesser,
      GGreaterDynamicsttp, GLesserDynamicsttp, ΣGreaterDynamicsttp, ΣLesserDynamicsttp)
    end
  end

  # Set retarded and advanced green function and self energy
  ΣRetardedDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 
  ΣAdavancedDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 

  GRetardedDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 
  GAdavancedDynamicsttp = zeros(ComplexF64, nt, nt, 2, 2) 

  for a in 1:2
    for b in 1:2
      ΣRetardedDynamicsttp[:,:,a,b] = LowerTriangular(ΣGreaterDynamicsttp[:,:,a,b]) - LowerTriangular(ΣLesserDynamicsttp[:,:,a,b])
      ΣAdavancedDynamicsttp[:,:,a,b] = UpperTriangular(ΣLesserDynamicsttp[:,:,a,b]) - UpperTriangular(ΣGreaterDynamicsttp[:,:,a,b])

      GRetardedDynamicsttp[:,:,a,b] = LowerTriangular(GGreaterDynamicsttp[:,:,a,b]) - LowerTriangular(GLesserDynamicsttp[:,:,a,b])
      GAdavancedDynamicsttp[:,:,a,b] = UpperTriangular(GLesserDynamicsttp[:,:,a,b]) - UpperTriangular(GGreaterDynamicsttp[:,:,a,b])
      for i in 1:nt
        ΣRetardedDynamicsttp[i,i,a,b]  *= 0.5
        ΣAdavancedDynamicsttp[i,i,a,b] *= 0.5

        GRetardedDynamicsttp[i,i,a,b]  *= 0.5
        GAdavancedDynamicsttp[i,i,a,b] *= 0.5
      end
    end
  end

  return J, nt, Δt, f ,GGreaterDynamicsttp   ,GLesserDynamicsttp    ,ΣGreaterDynamicsttp   ,ΣLesserDynamicsttp    ,GRetardedDynamicsttp  ,GAdavancedDynamicsttp ,ΣRetardedDynamicsttp  ,ΣAdavancedDynamicsttp
  

end

function thermal_spectrum_F(
  ΣKDynamicsttp::Array{ComplexF64, 4},
  GRetardedDynamicsttp::Array{ComplexF64, 4}, ΣRetardedDynamicsttp::Array{ComplexF64, 4},
  β::Float64, Λω::Float64, Λt::Float64, nt::Int64, k_index::Int64)
   
     

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
  a2=-Λt+0.00000001;
  b2= Λt;
  # nt = 2001;
  equaldTabscissae = collect(0:1/(nt-1):1)
  Δt=Λt/((nt-1)/2);
  equaldTabscissaet=[(b2-a2)*equaldTabscissae[i]+a2 for i in 1:nt];
  # equaldTweightst=Λt*equaldTweights/sum(equaldTweights);

  abscissaeω = vcat((abscissaeω[(1 + nω)÷2:nω] .- Λω), [0.00000000001], (Λω .+ abscissaeω[1 : (1 + nω)÷2]))
  weightsω =   vcat((weightsω[(1 + nω)÷2:nω]), [0.], (weightsω[1 : (1 + nω)÷2]))

  weightsω = 2*Λω * (weightsω/ sum(weightsω))

  weightst= ones(nt - 2);
  weightst =   vcat([0.5], weightst, [0.5])
  weightst = 2*Λt * (weightst/ sum(weightst))

  invEqualdTTrans= Array{ComplexF64, 2}(undef, length(equaldTabscissaet), length(abscissaeω))
  for i in 1:length(equaldTabscissaet)
      for j in 1:length(abscissaeω)
        invEqualdTTrans[i, j] = weightsω[j] * exp(- 1.0im*equaldTabscissaet[i]*abscissaeω[j])
      end
  end


  EqualdTTrans= Array{ComplexF64, 2}(undef, length(abscissaeω), length(equaldTabscissaet))
  for i in 1:length(abscissaeω)
      for j in 1:length(equaldTabscissaet)
        EqualdTTrans[i, j] = weightst[j] * exp(1.0im*abscissaeω[i]*equaldTabscissaet[j])
      end
  end

  nF = [1/(exp(β*ω) + 1) for ω in abscissaeω]
  reverseNF = [1/(exp(-β*ω) + 1) for ω in abscissaeω]

  cutoff_length = (nt - 1) ÷ 2

  GRetardedthermal = Array{ComplexF64, 3}(undef, nt, 2, 2)
  GRomegathermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)
  ρωthermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)

  ΣKthermal = Array{ComplexF64, 3}(undef, nt, 2, 2)
  ΣKomegathermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)
  ΣRmΣAthermal = Array{ComplexF64, 3}(undef, nt, 2, 2)
  ΣRmΣAomegathermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)
  Fomegathermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)

  # TODO pay attention to the order of GGreaterthermal(t)
  # Wrong order, FIXED
  for a in 1:2
    for b in 1:2
      GRetardedthermal[:,a,b] = vcat(
        GRetardedDynamicsttp[k_index - cutoff_length:k_index - 1,k_index,a,b],
        [GRetardedDynamicsttp[k_index, k_index,a,b]],
        reverse(GRetardedDynamicsttp[k_index, k_index - cutoff_length:k_index - 1,a,b]) 
          )
    end
  end

  for a in 1:2
      for b in 1:2
        GRomegathermal[:,a,b] = ((EqualdTTrans*GRetardedthermal[:,a,b]))
        # ρωthermal[1:500,a,b] = zeros(ComplexF64, 500)
      end
  end

  ####### WARNING: pay attention to the 2π!!!! the fourier transform doesn't include 2π
  for i in 1:length(abscissaeω)
    ρωthermal[i,:,:] = (GRomegathermal[i,:,:] - adjoint(GRomegathermal[i,:,:]))/(-2π*1.0im)
  end

  # Distribution function calculation
  for a in 1:2
    for b in 1:2
      ΣKthermal[:,a,b] = vcat(
        ΣKDynamicsttp[k_index - cutoff_length:k_index - 1,k_index,a,b],
        [ΣKDynamicsttp[k_index, k_index,a,b]],
        reverse(ΣKDynamicsttp[k_index, k_index - cutoff_length:k_index - 1,a,b]) 
          )
    end
  end

  for a in 1:2
    for b in 1:2
      ΣKomegathermal[:,a,b] = ((EqualdTTrans*ΣKthermal[:,a,b]))
    end
  end

  for a in 1:2
    for b in 1:2
      ΣRmΣAthermal[:,a,b] = vcat(
        ΣRetardedDynamicsttp[k_index - cutoff_length:k_index - 1,k_index,a,b],
        [ΣRetardedDynamicsttp[k_index, k_index,a,b]],
        reverse(ΣRetardedDynamicsttp[k_index, k_index - cutoff_length:k_index - 1,a,b]) 
          )
    end
  end

  for a in 1:2
    for b in 1:2
      ΣRmΣAomegathermal[:,a,b] = ((EqualdTTrans*ΣRmΣAthermal[:,a,b]))
    end
  end

  eta = 10^(-5)
  for i in 1:length(abscissaeω)
    ΣRmΣAomegathermal[i,:,:] = (ΣRmΣAomegathermal[i,:,:] - adjoint(ΣRmΣAomegathermal[i,:,:])) .+ eta*1im
  end
  Fomegathermal = ΣKomegathermal ./ ΣRmΣAomegathermal

  # ΣRetardedDynamicsttp - ΣRetardedDynamicsttp

  # temp_plot = Plots.plot(abscissaeω, hcat(real(ρωthermal[:,1,1]), imag(ρωthermal[:,1,1])))
  # savefig(temp_plot, "temp_plot.pdf")

  # JLD.save("ρω_thermal_kindex$(k_index).jld","arr", ρωthermal, "omega", abscissaeω, "k_index_plot", k_index)
  return ρωthermal, Fomegathermal

end


function thermal_Gt(GLesserDynamicsttp::Array{ComplexF64, 4}, β::Float64, Λω::Float64, Λt::Float64, nt::Int64, k_index::Int64)


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
  a2=-Λt+0.00000001;
  b2= Λt;
  # nt = 2001;
  equaldTabscissae = collect(0:1/(nt-1):1)
  Δt=Λt/((nt-1)/2);
  equaldTabscissaet=[(b2-a2)*equaldTabscissae[i]+a2 for i in 1:nt];

  cutoff_length = (nt - 1) ÷ 2
  # equaldTweightst=Λt*equaldTweights/sum(equaldTweights);

  abscissaeω = vcat((abscissaeω[(1 + nω)÷2:nω] .- Λω), [0.00000000001], (Λω .+ abscissaeω[1 : (1 + nω)÷2]))
  weightsω =   vcat((weightsω[(1 + nω)÷2:nω]), [0.], (weightsω[1 : (1 + nω)÷2]))

  weightsω = 2*Λω * (weightsω/ sum(weightsω))

  weightst= ones(nt - 2);
  weightst =   vcat([0.5], weightst, [0.5])
  weightst = 2*Λt * (weightst/ sum(weightst))

  invEqualdTTrans= Array{ComplexF64, 2}(undef, length(equaldTabscissaet), length(abscissaeω))
  for i in 1:length(equaldTabscissaet)
      for j in 1:length(abscissaeω)
        invEqualdTTrans[i, j] = weightsω[j] * exp(- 1.0im*equaldTabscissaet[i]*abscissaeω[j])
      end
  end


  EqualdTTrans= Array{ComplexF64, 2}(undef, length(abscissaeω), length(equaldTabscissaet))
  for i in 1:length(abscissaeω)
      for j in 1:length(equaldTabscissaet)
        EqualdTTrans[i, j] = weightst[j] * exp(1.0im*abscissaeω[i]*equaldTabscissaet[j])
      end
  end

  nF = [1/(exp(β*ω) + 1) for ω in abscissaeω]
  reverseNF = [1/(exp(-β*ω) + 1) for ω in abscissaeω]


  Glessomegathermal = Array{ComplexF64, 3}(undef, length(abscissaeω), 2, 2)

  # Plot thermal Gt
  Glessthermal = Array{ComplexF64, 3}(undef, nt, 2, 2)

  # TODO pay attention to the order of GGreaterthermal(t)
  # Wrong order, FIXED
  for a in 1:2
    for b in 1:2
      Glessthermal[:,a,b] = vcat(
        GLesserDynamicsttp[k_index - cutoff_length:k_index - 1,k_index,a,b],
        [GLesserDynamicsttp[k_index, k_index,a,b]],
        reverse(GLesserDynamicsttp[k_index, k_index - cutoff_length:k_index - 1,a,b]) 
          )
    end
  end

  for a in 1:2
      for b in 1:2
        Glessomegathermal[:,a,b] = ((EqualdTTrans*Glessthermal[:,a,b]))
        # ρωthermal[1:500,a,b] = zeros(ComplexF64, 500)
      end
  end

  return Glessthermal, Glessomegathermal

end