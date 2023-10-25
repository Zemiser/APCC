function DelayAPCC = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i)

TEnd_r = zeros(r,1); % the end delay of each set

% update the delay of each subtask for each worker
for index4 = 2:r
    VectorTm(index4,:) = VectorTm(index4,:) + VectorTm(index4 - 1,:);
end
[SortVectorT,~] = sort(VectorTm,2);
for index5 = 1:r % to obtain the delay of the complete task
    TEnd_r(index5) = SortVectorT(index5,VectorEta_i(index5));
end

DelayAPCC = max(TEnd_r); % to obtain the delay of the complete task

end