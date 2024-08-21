module SpikeSeqProcessor
    function EnumerateSpikedNeuron( t1, t2, SpikeTime::Vector, SpikeNeuron::Vector, N )
        i=1
        spike_binary = zeros( UInt32, N)
        for i=1:length(SpikeTime)
            if t1 <= SpikeTime[i] && SpikeTime[i] <= t2
                spike_binary[UInt(SpikeNeuron[i])] == 1
            end
        end
        return spike_binary
    end
    
    function EstimateX( SpikeTime::Vector, SpikeNeuron::Vector, ϵ1, ϵ2, T, N )
        num_bin = UInt32(floor((T-ϵ2)/ϵ1))
        numEvent_Yi1_Yj1 = zeros( UInt32, N, N)
        numEvent_Yi0_Yj1 = zeros( UInt32, N, N)
        numEvent_Yi1_Yj0 = zeros( UInt32, N, N)
        numEvent_Yi0_Yj0 = zeros( UInt32, N, N)
        for bin=1:num_bin
            spike_binary1 = EnumerateSpikedNeuron( ϵ1*(bin-1), ϵ1*bin   , SpikeTime, SpikeNeuron, N )
            spike_binary2 = EnumerateSpikedNeuron( ϵ1*bin    , ϵ1*bin+ϵ2, SpikeTime, SpikeNeuron, N )
            for i=1:N
                for j=1:N
                    if     spike_binary1[j] == 1 && spike_binary2[i] == 1
                        numEvent_Yi1_Yj1[i, j] += 1
                    elseif spike_binary1[j] == 0 && spike_binary2[i] == 1
                        numEvent_Yi1_Yj0[i, j] += 1
                    elseif spike_binary1[j] == 1 && spike_binary2[i] == 0
                        numEvent_Yi0_Yj1[i, j] += 1
                    elseif spike_binary1[j] == 0 && spike_binary2[i] == 0
                        numEvent_Yi0_Yj0[i, j] += 1
                    end
                end
            end
        end
        # 要素ごとの割り算
        # DevideError を例外処理実装する
        X = numEvent_Yi1_Yj1 ./ (numEvent_Yi1_Yj1 .+ numEvent_Yi1_Yj0) .- numEvent_Yi1_Yj0 ./ (numEvent_Yi1_Yj0 .+ numEvent_Yi0_Yj0)
        for i=1:N
            X[i, i] = 0.
        end
        return X
    end
    
    
    function Extract_SubPopulation( SpikeTime::Vector, SpikeNeuron::Vector, sub_population::Vector )
        SpikeTime_SubPop = []
        SpikeNeuron_SubPop = []
        for i=1:length(SpikeTime)
            if sub_population[SpikeNeuron[i]] != 0
                push!(SpikeTime_SubPop, SpikeTime[i])
                push!(SpikeNeuron_SubPop, sub_population[SpikeNeuron[i]])
            end
        end
        return SpikeTime_SubPop, SpikeNeuron_SubPop
    end
    
    
    function Extract_SubSequence( SpikeTime::Vector, SpikeNeuron::Vector, t_start, t_end )
        i = 1
        SpikeTime_SubSeq = []
        SpikeNeuron_SubSeq = []
        while SpikeTime[i] < t_start 
            i += 1
        end
        while t_start <= SpikeTime[i] && SpikeTime[i] <= t_end
            push!(SpikeTime_SubSeq, SpikeTime[i])
            push!(SpikeNeuron_SubSeq, SpikeNeuron[i])
            i += 1
        end
        return SpikeTime_SubSeq, SpikeNeuron_SubSeq
    end
end
