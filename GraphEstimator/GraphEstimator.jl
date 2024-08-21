using LinearAlgebra
using Random
include("SpikeSequence_Processor.jl")

function calc_α( X, Λ, N_recorded)
    temp = inv( I(N_recorded)-Λ) - I(N_recorded)
    top = 0.
    bottom = 0.
    for i=1:N_recorded
        for j=1:N_recorded
            if i != j
                top += X[i, j] * temp[i, j]
                bottom += temp[i,j]^2
            end
        end
    end
    return top / bottom
end

function calc_D( α, Λ, N_recorded )
    D = α.* ( inv(I(N_recorded)-Λ) - I(N_recorded))
    # take only diagonal elements
    for i=1:N_recorded
        for j=1:N_recorded
            if i != j
                D[i,j] = 0.
            end
        end
    end
    return D
end


function calc_Λ( α, X, D , M, N_recorded )
    Λ = (-1) .* inv(( 1/α).*(X + D) + I(N_recorded) )
    for i=1:N_recorded
        Λ[i, i] = 0.
    end

    #reduce element following the condtion for sparsity 
    Λ_sorted= sort(vec(Λ), by=abs, rev=true)
    ζ= abs(Λ_sorted[M])
    for i=1:N_recorded
        for j=1:N_recorded
            if abs(Λ[i, j]) < ζ
                Λ[i, j] = 0.
            end
        end
    end
    return Λ
end


function minimizeJ( X, r, M, num_iteration, N_recorded )
    # initial state
    α = rand()
    D  = r*I(N_recorded)
    Λ = zeros(N_recorded, N_recorded)

    for i=1:num_iteration
        Λ = calc_Λ( α, X, D , M, N_recorded)
        @printf("end Λ:%d\n", i)
        #println(Λ)
        α = calc_α( X , Λ, N_recorded )
        @printf("end α:%d\n", i)
        #println(α)
        D  = calc_D( α, Λ, N_recorded )
        @printf("end D:%d\n", i)
        #println(D)
    end

    return Λ
end


function StARS!( SpikeTime::Vector, SpikeNeuron::Vector, T, b, num_iteration, γ, K, N_recorded, ϵ1, ϵ2, r )
    SubSeq_start = (T-b) * rand(K)
    ξ_Mij = zeros(Float64, N_recorded, N_recorded)
    X_kij = zeros(Float64, K, N_recorded, N_recorded)

    # generate X for each subsequence
    for k=1:K
            SpikeTime_SubSeq, SpikeNeuron_SubSeq = SpikeSeqProcessor.Extract_SubSequence( SpikeTime, SpikeNeuron, SubSeq_start[k], SubSeq_start[k]+b )
            X_kij[k, :, :] = SpikeSeqProcessor.EstimateX( SpikeTime_SubSeq, SpikeNeuron_SubSeq, ϵ1, ϵ2, T, N_recorded)
    end


    for M=1:N_recorded*(N_recorded-1) # この繰返し範囲は，速度向上のために絞れる可能性がある
        ϕ_Mij = zeros(Float64, N_recorded, N_recorded)
        for k=1:K
            Λ = minimizeJ( X_kij[K, :, :], r, M, num_iteration, N_recorded )
            for i=1:N_recorded
                for j=1:N_recorded
                    if Λ[i, j] < -1e-10 &&  1e-10 < Λ[i, j] 
                        ϕ_Mij[i, j] += 1
                    end
                end
            end
        end

        # calculate probability and its variance
        ϕ_Mij = ϕ_Mij ./ K
        ξ_Mij = 2 * ϕ_Mij .* (1 .- ϕ_Mij)
        println(ξ_Mij)

        # calculate average of instability values
        D_M = 0.
        for i=1:N_recorded
            for j=1:N_recorded
                D_M += ξ_Mij[i, j]
            end
        end
        D_M = D_M / (N_recorded*(N_recorded-1))
        @printf("M:%d D_M:%.5f\n", M, D_M)

        # M is the maximum number which satisfy the condition that for all M' (<= M), Db(M') < γ
        if D_M > γ
            return M-1
        end
    end

    return N_recorded*(N_recorded-1) #最後までfor回してしまったら，どこかおかしい
end



Random.seed!(1)
# record timeseries
record_neuron = rand(N)
record_prob = 0.8
N_recorded = 0
for i=1:N
    if record_neuron[i] <= 0.8
        record_neuron[i] = 1
        global N_recorded += 1
    else
        record_neuron[i] = 0
    end
end

# rename neuron number i to new_name
# not recorded neurons are represented as 0 in rename_list
rename_list = zeros(N)
new_name = 1
for i=1:N
    if record_neuron[i] == 1
        push!(rename_list, new_name)
        rename_list[i] = new_name
        global new_name += 1
    end
end
println(record_neuron)
println(rename_list)


SpikeTime_SubPopulation, SpikeNeuron_SubPopulation = SpikeSeqProcessor.Extract_SubPopulation( SpikeTime, SpikeNeuron, rename_list )

## hyper parameters
ϵ1 = 3.0 #msec
ϵ2 = 3.0 #msec
r = 0.3
K = 100
γ= 0.2
b = 3000
num_iteration = 5
println("StARS started")
M = StARS!( SpikeTime_SubPopulation, SpikeNeuron_SubPopulation, T, b, num_iteration, γ, K, N_recorded, ϵ1, ϵ2, r)
println("StARS ended")

# optimization
X = SpikeSeqProcessor.EstimateX( SpikeTime_SubPopulation, SpikeNeuron_SubPopulation, ϵ1, ϵ2, T, N_recorded )
println("Optimization started")
Λ = minimizeJ( X, r, M, num_iteration, N_recorded )
println("Optimization ended")
Λ_ans = zeros(N_recorded, N_recorded)
for i=1:N
    for j=1:N
        if rename_list[i] != 0 && rename_list[j] != 0
            Λ_ans[Int(rename_list[i]), Int(rename_list[j])] = W[i, j]
        end
    end
end

p1 = plot(Λ, clim=(-1,1))
p2 = plot(Λ_ans, clim=(-1,1))
plot(p1, p2, layout=(2,1))
