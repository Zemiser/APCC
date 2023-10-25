function [VectorKStar_i, Flag, tStar] = T3BasedK_i(H,N,d,L,mu,alpha,r,Delta1)
% The result of K_i in this function can be non-integer

Flag = 1; %1: solvable 0: unsolvable
VectorKStar_i = zeros(1,r); % K_i = VectorKStar_i(i+1);
VectorEtaStar_i = zeros(1,r); % K_i = VectorKStar_i(i+1);
tStar = alpha * r + r / mu; %initialization
KPrimeMax = floor((N - d*(L-1)-1)/d)-1;

if H == r*d*(L+1) + r
    VectorKStar_i = 2*ones(1,r);
    tStar = alpha*r + r/mu/(N - (d*(L+1)+1) + 1);
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
    tStar = alpha*r + r/mu/(N - VectorEtaStar_i + 1);
    return
end

tStarTemp1 = r * alpha; % Minimum required delay
if  (r*(N+1) - H) <= 0 || SumForOptimalK_T3(tStarTemp1,alpha,r) <= (mu*r*(N+1) - mu*H) % To ensure solvability
    Flag = 0;
    return 
end

while SumForOptimalK_T3(tStarTemp1,alpha,r) > (mu*r*(N+1) - mu*H) %Traversal to determine the value range of tStar
    tStarTemp1 = tStarTemp1 * 2;
end
tStarTemp2 = tStarTemp1/2; %tStarTemp1 is the maximum of the interval and tStarTemp2 is the minimum
tStar = (tStarTemp2 + tStarTemp1)/2;
while tStarTemp1 - tStarTemp2 >= Delta1
    if SumForOptimalK_T3(tStar,alpha,r) >= (mu*r*(N+1) - mu*H)
        tStarTemp2 = tStar;
    else
        tStarTemp1 = tStar;
    end
    tStar = (tStarTemp2 + tStarTemp1)/2;
end

% To obtain the K_i
for index2 = 0:r-1
    VectorEtaStar_i(index2+1) = N + 1 - (index2 + 1)/mu/(tStar - alpha*(index2 + 1));
    VectorKStar_i(index2+1) = (VectorEtaStar_i(index2+1) - d*(L-1)-1)/d;
end
if VectorKStar_i(end) < 2 % ensure that K_i >= 2
    tStarTemp3 = alpha*r + r/mu/(N - (d*(L+1)+1) + 1);% the tStar corresponding to the last set which has K_(r-1)=2, note that tStarTemp3 is larger than alpha*r
    [VectorKStar_iTemp, Flag, tStarTemp4] = T3BasedK_i(H-d*(L+1)-1,N,d,L,mu,alpha,r-1,Delta1);
    tStar = max(tStarTemp3,tStarTemp4);
    if tStar <= alpha*r
        Flag = 0;
        return
    end
    VectorKStar_i = [VectorKStar_iTemp,2];
end
end
