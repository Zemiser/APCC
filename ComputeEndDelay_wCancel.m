function DelayAPCC = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i)

TEnd_r = zeros(r,1); % the end delay of each set

% update the delay of each subtask for each worker
for index4 = 1:r
    if index4 == 1
        [SortVectorTmTemp,SortIndex1] = sort(VectorTm(index4,:));
        TEnd_r(index4) = SortVectorTmTemp(VectorEta_i(index4)); % the end delay of this set
        VectorTm(index4,SortIndex1(VectorEta_i(index4)+1:end)) = TEnd_r(index4); % update the end delay of this set for other workers
    else
        VectorTm(index4,:) = VectorTm(index4,:) + VectorTm(index4 - 1,:);
        [SortVectorTmTemp,SortIndex2] = sort(VectorTm(index4,:));
        TEnd_r(index4) = SortVectorTmTemp(VectorEta_i(index4));% the end delay of this set
        for index5 = SortIndex2(VectorEta_i(index4)+1:end) % update the end delay of this set for other workers
            VectorTm(index4,index5) = max([VectorTm(index4 - 1,index5),TEnd_r(index4)]);
        end
    end
end

DelayAPCC = max(TEnd_r);% to obtain the delay of the complete task

end