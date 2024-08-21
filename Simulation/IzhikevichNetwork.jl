using Base: @kwdef
using Parameters: @unpack
using Plots
using Printf
using Random

@kwdef struct IZParameter{FT}
    a::FT = 0.1
    b::FT = 0.25
    c::FT = -65
    d::FT = 8
    vpeak::FT = 30
    noise_freq::FT = 10.0 #Hz
    noise_strength::FT = 10.
    exsyn_strength::FT = 10.
    insyn_strength::FT = 10.
    #exsyn_strength::FT = 0.
    #insyn_strength::FT = 0.
    τ_noise::FT = 1.
    τ_exsyn::FT = 2.
    τ_insyn::FT = 2.

end

@kwdef mutable struct IZ{FT}
    param::IZParameter = IZParameter{FT}()
    N_IN::UInt32
    N_EX::UInt32
    N::UInt32
    v::Vector{FT} = fill(param.c, N)
    u::Vector{FT} = zeros(N)
    I_insyn::Vector{FT} = zeros(N)
    I_exsyn::Vector{FT} = zeros(N)
    spike::Vector{Bool} = zeros(Bool, N)
end

function update!( variable::IZ, param::IZParameter, I_noise::Vector, dt, spike_time::Vector, spike_neuron::Vector, time )
    @unpack N, v, u, I_insyn, I_exsyn, spike = variable
    @unpack a, b, c, d, vpeak, noise_freq, noise_strength, exsyn_strength, insyn_strength = param

    @inbounds for i=1:N
    end

    @inbounds for i=1:N
        v[i] += dt*(0.04*v[i]^2+5*v[i]+140-u[i]+I_insyn[i]+I_exsyn[i]+I_noise[i])
        u[i] += dt*(a*(b*v[i]-u[i]))
    end

    @inbounds for i=1:N
        spike[i] = v[i] >= vpeak
        if spike[i] == true
            push!(spike_time, time)
            push!(spike_neuron, i)
        end
        v[i] = ifelse(spike[i], c, v[i])
        u[i] += ifelse(spike[i], d, 0)
    end
end

function calculate_synaptic_current( variable::IZ, param::IZParameter, W::Matrix )
    @unpack N, N_EX, N_IN, I_insyn, I_exsyn, spike = variable
    @unpack noise_freq, noise_strength, exsyn_strength, insyn_strength, τ_noise, τ_insyn, τ_exsyn = param

    for i=1:N
        I_exsyn[i] = exp(-dt/τ_exsyn)I_exsyn[i]
        I_insyn[i] = exp(-dt/τ_insyn)I_insyn[i]
        I_noise[i] = exp(-dt/τ_noise)I_noise[i]
    end


    for j=1:N_EX
        if spike[j] == true
            for i=1:N
                I_exsyn[i] += exsyn_strength*W[i, j]
            end
        end
    end

    for j=N_EX+1:N_IN
        if spike[j] == true
            for i=1:N
                I_insyn[i] += insyn_strength*W[i, j]
            end
        end
    end

    for i=1:N
        if rand() < noise_freq * 0.001 * dt
            I_noise[i] += noise_strength
        end
    end
end



Random.seed!(1)
T = 50000
dt = 0.01f0
nt = UInt32(T/dt)
N_EX = 40
N_IN = 10
N = N_EX + N_IN
t = Array{Float32}(1:nt)*dt
neurons = IZ{Float32}(N=N, N_IN=N_IN, N_EX=N_EX)

# determine synaptic connections
W = rand(N, N)
connection = rand(N, N)
connection_prob = 0.1
for j=1:N
    W[j, j] = 0
    for i=1:N
        W[i, j] = (connection[i, j] < connection_prob) ? W[i, j] : 0
    end
end
for j=N_EX+1:N
    for i=1:N
        W[i, j] = (-1)*W[i, j]
    end
end

SpikeTime = []
SpikeNeuron = []

I_noise = zeros(Float32, N)
@time for i=1:nt
    update!(neurons, neurons.param, I_noise, dt, SpikeTime, SpikeNeuron, t[i])
    calculate_synaptic_current(neurons, neurons.param, W)
end

scatter(SpikeTime, SpikeNeuron, xlims=(0,1000), ylims=(0,50))
#savefig("fig1.png")
