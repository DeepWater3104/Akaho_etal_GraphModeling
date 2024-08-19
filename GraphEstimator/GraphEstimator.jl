using LinearAlgebra
using Random

function EnumerateSpikedNeuron( t1, t2, SpikeTime::Vector, SpikeNeuron::Vector, N )
    i=1
    spike_binary = zeros( UInt32, N)
    while SpikeTime[i] < t1
        i++
    end
    while t1 <= SpikeTime[i] && SpikeTime[i] <= t2
        spike_binary[SpikeNeuron[i]] == 1
        i++
    end
    return spike_binary
end

function EstimateX_from_SpikeSequence( spike_time::Vector, spike_neuron::Vector, ϵ1, ϵ2, T, N )
    num_bin = UInt32((T-ϵ2)/ϵ1)
    numEvent_Yi1_Yj1 = zeros( UInt32, N, N)
    numEvent_Yi0_Yj1 = zeros( UInt32, N, N)
    numEvent_Yi1_Yj0 = zeros( UInt32, N, N)
    numEvent_Yi0_Yj0 = zeros( UInt32, N, N)
    for bin=1:num_bin
        spike_binary1 = EnumerateSpikedNeuron( ϵ1*(bin-1), ϵ1*bin   , spike_time, spike_neuron )
        spike_binary2 = EnumerateSpikedNeuron( ϵ1*bin    , ϵ1*bin+ϵ2, spike_time, spike_neuron )
        for i=1:N
            for j=1:N
                if     spike_binary1[j] == 1 && spike_binary2[i] == 1
                    numEvent_Yi1_Yj1[i, j] ++
                elseif spike_binary1[j] == 0 && spike_binary2[i] == 1
                    numEvent_Yi1_Yj0[i, j] ++
                elseif spike_binary1[j] == 1 && spike_binary2[i] == 0
                    numEvent_Yi0_Yj1[i, j] ++
                elseif spike_binary1[j] == 0 && spike_binary2[i] == 0
                    numEvent_Yi0_Yj0[i, j] ++
                end
            end
        end
    end
    # 要素ごとの割り算
    X = numEvent_Yi1_Yj1 ./ (numEvent_Yi1_Yj1 .+ numEvent_Yi1_Yj0) .- numEvent_Yi1_Yj0 ./ (numEvent_Yi1_Yj0 .+ numEvent_Yi0_Yj0)
    return X
end


function Extract_SubPopulation( spike_time::Vector, spike_neuron::Vector, sub_population::Vector )

end


function Extract_SubSequence( spike_time::Vector, spike_neuron::Vector, t1, t2 )
    
end


function calc_α( X, Λ)
    temp = inv( I(size(X,1))-Λ) - I(size(X,1))
    for i=1:N
        for j=1:N
            if i != j
                top += X[i, j] * temp[i, j]
                bottom += X[i, j] * temp[i,j]^2
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


function calc_Λ_fixedM( α, X, D , M )
    Λ = (-1)* inv(( 1/α)*(X + D) + I(size(X,1)) )
    for i=1:N
        Λ[i, i] = 0.
    end
    #reduce element
    
    return Λ
end

function minimizeJ_fixedM( num_iteration )
    Λ = calc_Λ_fixedM( α, X, D )
    α = calc_α( X , Λ )
    D  = calc_D ( α, Λ )
    return Λ
end

function StARS!( SpikeTime::Vector, SpikeNeuron::Vector, T, num_iteration )
    # derive subsequences
    subseq = N*rand(2, K)
    ψMkij = zeros(N,K,N,N)
    for M=1:N
        for k=1:K
            Xsub = extract_subseq( SpikeTime, SpikeNeuron, subseq[1,k], subseq[2,k] )
            minimizeJ_fixedM( num_iteration )
            ψMkij [M, k, :, :] = 
        end
    end
end

function calc_Λ( α, X, D )
    Λ = (-1)* inv(( 1/α)*(X + D) + I(size(X,1)) )
    for i=1:N
        Λ[i, i] = 0.
    end
    StARS!( Λ)
    return Λ
end

function minimizeJ( X, num_iteration )
    # initial state
    α = rand()
    D  = r*I(N)
    Λ = zeros(N, N)

    for i=1:num_iteration
        Λ = calc_Λ( α, X, D )
        α = calc_α( X , Λ )
        D  = calc_D ( α, Λ )
    end
    return Λ
end


# hyper parameters
Random.seed!(1)
ϵ1 = 3.0 #msec
ϵ2 = 3.0 #msec
r = 1.0
K = 100
γ= 
num_iteration = 5
record_neuron = []

# optimization
SpikeTime_SubPopulation, SpikeTime_SubPopulation = extract_subPopulation( SpikeTime, SpikeNeuron, record_neuron )
X = EstimateX_from_SpikeSequence( SpikeTime, SpikeNeuron, ϵ1, ϵ2, T, N )
Λ = minimizeJ( X, num_iteration )
