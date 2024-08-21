using LinearAlgebra
using Random
include("SpikeSequence_Processor.jl")

function calc_α( X, Λ)
    temp = inv( I(size(X,1))-Λ) - I(size(X,1))
    top = 0.
    bottom = 0.
    for i=1:N
        for j=1:N
            if i != j
                top += X[i, j] * temp[i, j]
                bottom += temp[i,j]^2
            end
        end
    end
    return top / bottom
end

function calc_D( α, Λ)
    D = α* ( inv(I(size(Λ,1))-Λ) - I(size(Λ,1)) )
    # take only diagonal elements
    for i=1:size(D ,1)
        for j=1:size(D ,1)
            if i != j
                D[i,j] = 0.
            end
        end
    end
    return D
end


function calc_Λ( α, X, D , M )
    Λ = (-1)* inv(( 1/α)*(X + D) + I(size(X,1)) )
    for i=1:N
        Λ[i, i] = 0.
    end

    #reduce element following the condtion for sparsity 
    Λ_sorted= sort(vec(Λ), by=abs)
    for i=1:N
        for j=1:N
            if abs(Λ[i, j]) < abs(Λ_sorted[M])
                Λ[i, j] = 0.
            end
        end
    end
    return Λ
end


function minimizeJ( X, r, M, num_iteration )
    # initial state
    α = rand()
    D  = r*I(N)
    Λ = zeros(N, N)

    for i=1:num_iteration
        Λ = calc_Λ( α, X, D , M)
        α = calc_α( X , Λ )
        D  = calc_D( α, Λ )
    end

    return Λ
end


function StARS!( SpikeTime::Vector, SpikeNeuron::Vector, T, b, num_iteration, γ, K, N, ϵ1, ϵ2, r )
    SubSeq_start = (T-b) * rand(K)
    ξ_Mij = zeros(Float64, N, N)
    X_kij = zeros(Float64, K, N, N)

    # generate X for each subsequence
    for k=1:K
            SpikeTime_SubSeq, SpikeNeuron_SubSeq = SpikeSeqProcessor.Extract_SubSequence( SpikeTime, SpikeNeuron, SubSeq_start[k], SubSeq_start[k]+b )
            X_kij[k, :, :] = SpikeSeqProcessor.EstimateX( SpikeTime_SubSeq, SpikeNeuron_SubSeq, ϵ1, ϵ2, T, N)
    end


    for M=1000:N*(N-1) # この繰返し範囲は，速度向上のために絞れる可能性がある
        ϕ_Mij = zeros(Float64, N, N)
        for k=1:K
            Λ = minimizeJ( X_kij[K, :, :], r, M, num_iteration )
            for i=1:N
                for j=1:N
                    if Λ[i, j] < -1e-10 &&  1e-10 < Λ[i, j] 
                        ϕ_Mij[i, j] += 1
                    end
                end
            end
        end

        # calculate probability and its variance
        ϕ_Mij = ϕ_Mij ./ K
        ξ_Mij = 2 * ϕ_Mij .* (1 .- ϕ_Mij)

        # calculate average of instability values
        D_M = 0.
        for i=1:N
            for j=1:N
                D_M += ξ_Mij[i, j]
            end
        end
        D_M = D_M / (N*(N-1))
        @printf("M:%d D_M:%.5f\n", M, D_M)

        # M is the maximum number which satisfy the condition that for all M' (<= M), Db(M') < γ
        if D_M > γ
            return M-1
        end
    end

    return N*(N-1) # 絶対なんかおかしい
end



function main()
    # record timeseries
    record_neuron = rand(N)
    record_prob = 0.8
    for i=1:N
        if record_neuron[i] <= 0.8
            record_neuron[i] = 1
        else
            record_neuron[i] = 0
        end
    end
    SpikeTime_SubPopulation, SpikeNeuron_SubPopulation = SpikeSeqProcessor.Extract_SubPopulation( SpikeTime, SpikeNeuron, record_neuron )

    # hyper parameters
    Random.seed!(1)
    ϵ1 = 3.0 #msec
    ϵ2 = 3.0 #msec
    r = 1.0
    K = 100
    γ= 0.2
    b = 10
    num_iteration = 5
    println("StARS started")
    M = StARS!( SpikeTime_SubPopulation, SpikeNeuron_SubPopulation, T, b, num_iteration, γ, K, N, ϵ1, ϵ2, r)
    println("StARS ended")

    # optimization
    X = SpikeSeqProcessor.EstimateX( SpikeTime_SubPopulation, SpikeNeuron_SubPopulation, ϵ1, ϵ2, T, N )
    println("Optimization started")
    Λ = minimizeJ( X, r, num_iteration )
    println("Optimization ended")
    # re-rename neuron number here
end

if contains( @__FILE__, PROGRAM_FILE )
    @time main()
end
