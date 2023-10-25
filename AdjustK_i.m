function VectorK_i = AdjustK_i(VectorK_iTemp,KPrimeMax)

r = length(VectorK_iTemp);
Delta = VectorK_iTemp(1) - (KPrimeMax + 1);
VectorK_iTemp(1) = KPrimeMax + 1;
for index1 = 2:r
    while (Delta > 0 && VectorK_iTemp(index1) < (KPrimeMax + 1))
        VectorK_iTemp(index1) = VectorK_iTemp(index1) + 1;
        Delta = Delta - 1;
    end
    if Delta == 0
        break;
    end
end

VectorK_i = VectorK_iTemp;
end
