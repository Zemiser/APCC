function Sum = SumForOptimalK_T2(tStar,mu,alpha,r)
Sum = 0;
for i = 0:r-1
    Sum = Sum + exp(-mu*(tStar/(i + 1) - alpha));
end
end
