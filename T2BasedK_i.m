function [VectorKStar_i, Flag, tStar] = T2BasedK_i(H,N,d,L,mu,alpha,r,Delta1)
% The result of K_i in this function can be non-integer

Flag = 1; %1: solvable 0: unsolvable
VectorKStar_i = zeros(1,r); % K_i = VectorKStar_i(i+1);
VectorEtaStar_i = zeros(1,r); % K_i = VectorKStar_i(i+1);
tStar = alpha * r + r / mu; %initialization
KPrimeMax = floor((N - d*(L-1)-1)/d)-1;

if H == r*d*(L+1) + r
    VectorKStar_i = 2*ones(1,r);
    tStar = alpha*r - r*log(1 - (d*(L+1) + 1)/N)/mu;
    return
end

% Iteration termination condition
if r == 1
    VectorEtaStar_i = H;
    VectorKStar_i = (VectorEtaStar_i - d*(L-1)-1)/d;
    if VectorKStar_i < 2 || VectorKStar_i > KPrimeMax + 1
%         disp('ERROR of T2BasedK_i!');
        Flag = 0;
        return
    end
    tStar = alpha*r - r*log(1 - VectorEtaStar_i/N)/mu;
    return
end

tStarTemp1 = alpha; % Minimum required delay 初始化的时候可以选取alpha而非r*alpha,
% 此时SumForOptimalK_T2(tStarTemp1,mu,alpha,r) 一定大于 (r - H/N)，
% 这时候求得的t可能不满足条件，但最后只会导致K_i不满足条件（K_i<2）,可以通过后续手段继续求解
if  (r - H/N) <= 0 || SumForOptimalK_T2(tStarTemp1,mu,alpha,r) <= (r - H/N) % To ensure solvability
    Flag = 0;
    return 
end

while SumForOptimalK_T2(tStarTemp1,mu,alpha,r) > (r - H/N) %Traversal to determine the value range of tStar
    tStarTemp1 = tStarTemp1 * 2;
end
tStarTemp2 = tStarTemp1/2; %tStarTemp1 is the maximum of the interval and tStarTemp2 is the minimum
tStar = (tStarTemp2 + tStarTemp1)/2;
while tStarTemp1 - tStarTemp2 >= Delta1
    if SumForOptimalK_T2(tStar,mu,alpha,r) >= (r - H/N)
        tStarTemp2 = tStar;
    else
        tStarTemp1 = tStar;
    end
    tStar = (tStarTemp2 + tStarTemp1)/2;
end

% To obtain the K_i
for index2 = 0:r-1
    VectorEtaStar_i(index2+1) = N*(1 - exp(-mu*(tStar/(index2 + 1) - alpha)));
    VectorKStar_i(index2+1) = (VectorEtaStar_i(index2+1) - d*(L-1)-1)/d;
end
if VectorKStar_i(end) < 2 % ensure that K_i >= 2
    tStarTemp3 = alpha*r - r*log(1 - (d*(L+1) + 1)/N)/mu;% the tStar corresponding to the last set which has K_(r-1)=2, note that tStarTemp3 is larger than alpha*r
    [VectorKStar_iTemp, Flag, tStarTemp4] = T2BasedK_i(H-d*(L+1)-1,N,d,L,mu,alpha,r-1,Delta1);
    tStar = max(tStarTemp3,tStarTemp4);
    if tStar <= alpha*r
        Flag = 0;
        return
    end
    VectorKStar_i = [VectorKStar_iTemp,2];
end
end
