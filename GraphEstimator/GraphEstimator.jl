using LinearAlgebra
using Random

function EnumerateSpikedNeuron( t1, t2, spike_time::Vector, spike_neuron::Vector, N )
    i=1
    spike_binary = zeros( UInt32, N)
    while spike_time[i] < t1
        i++
    end
    while t1 <= spike_time[i] && spike_time[i] <= t2
        spike_binary[spike_neuron[i]] == 1
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
    X = numEvent_Yi1_Yj1 / (numEvent_Yi1_Yj1 + numEvent_Yi1_Yj0) - numEvent_Yi1_Yj0 / (numEvent_Yi1_Yj0 + numEvent_Yi0_Yj0)
    return X
end


function StARS!()
end

function calc_Λ( α, X, D )
    Λ = (-1)* inv(( 1/α)*(X + D) + I(size(X,1)) )
    for i=1:N
        Λ[i, i] = 0.
    end
    StARS!( Λ)
    return Λ
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
    for i=1:size(Λ,1)
        for i=1:size(Λ,1)
            if i != j
                D[i,j] = 0.
            end
        end
    end
    return D
end


Random.seed!(1)
ϵ1 = 3.0 #msec
ϵ2 = 3.0 #msec
r = 1.0
num_iteration = 5
X = EstimateX_from_SpikeSequence( spike_time, spike_neuron, ϵ1, ϵ2, T, N )

# initial state
α = rand()
D  = r*I(N)
Λ = zeros(N, N)

# optimization
for i=1:num_iteration
    Λ = calc_Λ( α, X, D )
    α = calc_α( X , Λ )
    D  = calc_D ( α, Λ )
end
