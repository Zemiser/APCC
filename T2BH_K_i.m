function [VectorKStar_i, Flag, tStar] = T2BH_K_i(K_remaining,N,d,L,mu,alpha,r,tStarUB,VectorK_i_Rounding,index_r,LastK_i,Delta1,VectorKStar_i_T2)
% K_remaining: the remaining K in the search process
% index_r: the current index of searching set
% The results of K_i in this function are integers

% if r == 4 && K_remaining == 148 && index_r == r
%     disp('test point');
% end


Flag = 0; %1: solvable 0: unsolvable 2: equal to T2Based
VectorKStar_i = VectorK_i_Rounding; % K_i = VectorKStar_i(i+1);
% VectorEtaStar_i = zeros(1,index_r); % K_i = VectorKStar_i(i+1);
tStar = tStarUB; %initialization, to search the lower bound

KPrimeMax = floor((N - d*(L-1)-1)/d)-1;
% Iteration termination condition
if index_r == 1
    VectorKStar_i = K_remaining;
    VectorEtaStar_i = d*(VectorKStar_i + L -1) +1; 
    if (VectorKStar_i < LastK_i) || (VectorKStar_i > KPrimeMax + 1)
        Flag = 0;
        return
    end
    Flag = 1;
    tStar = alpha*index_r - index_r*log(1 - VectorEtaStar_i/N)/mu;
    return
end


minimumK_i = max(floor(VectorKStar_i_T2(index_r)) - 1,LastK_i);
maximumK_i = min(ceil(VectorKStar_i_T2(index_r)) + 1, KPrimeMax + 1);
countTimes = -1; % to count the number of optimal solutions
for K_i = minimumK_i:maximumK_i
    Eta_i = d*(K_i + L - 1) + 1;
    tStarTemp1 = alpha*index_r - index_r*log(1 - Eta_i/N)/mu;
    if tStarTemp1 > tStarUB
        continue;
    else
        H_remaining = d*(K_remaining - K_i + (index_r-1)*L - (index_r-1)) + (index_r-1);
        [VectorKStar_iTemp1, FlagTemp, tStarTemp2] = T2BasedK_i(H_remaining,N,d,L,mu,alpha,index_r-1,Delta1);
    end

    if (FlagTemp == 0) || (tStarTemp2 > tStarUB)
%     if tStarTemp2 > tStarUB
        continue;
    else
        [VectorKStar_iTemp2, Flag, tStarTemp3] = T2BH_K_i(K_remaining - K_i,N,d,L,mu,alpha,r,tStarUB,VectorK_i_Rounding,index_r-1,K_i,Delta1,VectorKStar_i_T2);
        VectorKStar_iTemp = [VectorKStar_iTemp2, K_i];
        tStarTemp = max(tStarTemp1,tStarTemp3);
    end
    
    if (Flag == 1) && (tStarTemp <= tStar)
        if tStarTemp < tStar || countTimes == -1
            tStar = tStarTemp;
            VectorKStar_i = VectorKStar_iTemp;
            countTimes = 0; % reset
        else
            countTimes = countTimes + 1;
        end
    end
end

if countTimes >= 0
    Flag = 1; % to show that there is at least one solution
end
if countTimes >= 1
    disp('More than one optimal solution!');
end

end