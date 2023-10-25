function VectorK_i = RoundingK_i(K,r,VectorK_iTemp,KPrimeMax)

VectorK_iTemp(VectorK_iTemp>=(KPrimeMax+1)) = (KPrimeMax+1); %可能存在rounding之后刚好比K_imax要大1的情况。
DeltaK = K - sum(VectorK_iTemp); % To ensure the sum of K_i equal to K, basic idea: keep K_i >= K_{i+1}
VectorK_i = VectorK_iTemp;
while DeltaK ~= 0
    for index2 = 1:r
        if DeltaK <= 0
            break;
        end
        if VectorK_i(index2) < KPrimeMax + 1
            VectorK_i(index2) = VectorK_i(index2) + 1;
            DeltaK = DeltaK - 1;
        end
    end
    for index2 = 1:r
        if DeltaK >= 0
            break;
        end
        if VectorK_i(r+1-index2) > 2
            VectorK_i(r+1-index2) = VectorK_i(r+1-index2) - 1;
            DeltaK = DeltaK + 1;
        end
    end
end

end