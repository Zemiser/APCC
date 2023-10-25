function [K_simu_opt_vec, t_simu_opt] = Brute_force_wo_Cancel(K_in_vec,r,K,N,d,L,...
    mu0,alpha0,Delta,K_i_max,RepeatTimes)

% (minK_i,maxK_i,VectorK_i,index,K_remaining,DelayAPCCMinTemp,VectorK_iMinTemp,RepeatTimes,N,d,L,K,mu0,alpha0,VectorK_i_T2BH,Delta2)
% search around the value of VectorK_i_T2BH

%% 获得所有满足要求的K_vec: K_search_mat
K_i_range_mat = cell(1,r);%表示K_vec的搜索空间
for index_i = 1:r
    range_min = max(K_in_vec(index_i) - Delta, 2);%搜索空间中的K_i要大于等于2
    range_max = min(K_in_vec(index_i) + Delta, K_i_max + 1); % 搜索空间中的K_i要小于等于 (K_i_max + 1)
    K_i_range_mat{index_i} = (range_min:range_max)';
end

range_len_vec = zeros(1,r); %每个K_i的取值范围的长度
for index_i = 1:r
    range_len_vec(index_i) = length(K_i_range_mat{index_i});
end

Num_temp_K_vec = 1; % 表示搜索空间中大概一共有多少个K_vec
for index_i = 1:r
    Num_temp_K_vec = Num_temp_K_vec*range_len_vec(index_i);
end

K_searchtemp_mat = zeros(Num_temp_K_vec,r); % 搜索空间中的所有K_vec
K_search_mat = zeros(1,r); %经过筛选后符合要求的K_vec
Num_K_vec = 1; %经过筛选后符合要求的K_vec的数量

% Progress_ratio = 0;
for index1 = 1:Num_temp_K_vec %当r很大时，这个范围会指数速度扩大

    %     if index1/Num_temp_K_vec - Progress_ratio >= 0.01
    %         Progress_ratio = index1/Num_temp_K_vec;
    %         disp(['M =',num2str(M),',(BF) Searching K Progress:',num2str(100*Progress_ratio),'%']);
    %     end

    Flag = 1; %表示这个K_searchtemp_mat(index1,:)是否符合要求
    K_i_index_vec = zeros(1,r); %表示K_search_mat(index1,:)这个K_vec中每个K_i在其相应取值范围中的index
    dividend = index1; %表示在逐步获得上一行代码index过程中的被除数，初始化为index1
    for index_i = 1:r
        divisor = 1; %divisor表示除数，即之后K_i取值范围长度之积
        if index_i < r %如果是最后一个K_r,divisor = 1；
            for index2 = (index_i+1):r
                divisor = divisor*range_len_vec(index2);
            end
        end
        multiplier = ceil(dividend/divisor); %倍数
        K_i_index_vec(index_i) = multiplier;
        if multiplier < 1
            disp('Error_BF');
            pause
        end
        dividend = dividend - (multiplier - 1)*divisor;
        K_searchtemp_mat(index1,index_i) = K_i_range_mat{index_i}(K_i_index_vec(index_i));

        %检验是否满足K_i从大到小变化
        if index_i > 1 && K_searchtemp_mat(index1,index_i - 1) < K_searchtemp_mat(index1,index_i)
            Flag = 0;
            break;
        end
    end

    if Flag == 0
        continue;
    end

    if K - sum(K_searchtemp_mat(index1,:)) ~= 0
        continue;
    end

    K_search_mat(Num_K_vec,:) = K_searchtemp_mat(index1,:); %这些是符合要求的K_vec
    Num_K_vec = Num_K_vec + 1;
end
Num_K_vec = Num_K_vec - 1; %最后减去多余的一个
K_search_mat = K_search_mat(1:Num_K_vec,:); %把多余的0删去


t_simu_opt = 10*(alpha0 + 1/mu0); %initialization
K_simu_opt_vec = K_in_vec;

%遍历所有K_in_vec附近符合要求的K_vec
for index1 = 1:Num_K_vec
    VectorEta_i = d*(K_search_mat(index1,:) + L -1) +1;
    SumDelayAPCC = 0;
    for indexRepeat = 1: RepeatTimes
        u1 = rand(r,N);
        VectorT = alpha0 / K - log(1-u1) / (mu0 * K); % the delay of each subtask for each worker
        DelayAPCCTemp = ComputeEndDelay_woCancel(VectorT,r,VectorEta_i);
        SumDelayAPCC = SumDelayAPCC + DelayAPCCTemp;
    end
    DelayAPCC = SumDelayAPCC / RepeatTimes;

    if DelayAPCC < t_simu_opt
        t_simu_opt = DelayAPCC;
        K_simu_opt_vec = K_search_mat(index1,:);
    end
end
end

