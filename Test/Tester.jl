function judge_zero(x)
    if -1e-10 < x && x < 1e-10
        return true
    else 
        return false
    end
end

# True Positive
TP = 0
# True Negative
TN = 0
# False Positive
FP = 0
# False Negative
FN = 0
TrueOrFalse = zeros(N_recorded, N_recorded)
for j=1:N_recorded
    for i=1:N_recorded
       if !judge_zero(Λ[i, j]) && !judge_zero(Λ_ans[i, j])
           global TP += 1
           TrueOrFalse[i, j] = true
       elseif judge_zero(Λ[i, j]) && judge_zero(Λ_ans[i, j])
           global TN += 1
           TrueOrFalse[i, j] = true
       elseif !judge_zero(Λ[i, j]) && judge_zero(Λ_ans[i, j])
           global FP += 1
           TrueOrFalse[i, j] = false
       elseif judge_zero(Λ[i, j]) && !judge_zero(Λ_ans[i, j])
           global FN += 1
           TrueOrFalse[i, j] = false
       else
           println("missed")
       end

    end
end

Sensitivity = TP / (TP + FN)
FP_rate = FP / (FP + TN)

heatmap(TrueOrFalse)
savefig("TrueOrFalse.png")

@printf("TP:%f TN:%f FP:%f FN:%f\n", TP, TN, FP, FN)
@printf("Sensitivity:%f\n", Sensitivity)
@printf("FP rate:%f\n", FP_rate)
