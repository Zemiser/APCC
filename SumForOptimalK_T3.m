function Sum = SumForOptimalK_T3(zStar,alpha,r)
Sum = 0;
for i = 0:r-1
    Sum = Sum + (i + 1)/(zStar - alpha*(i + 1));
end
end
