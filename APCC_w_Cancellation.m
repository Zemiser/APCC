clc;clear;close all;
%% For performance verification of APCC with cancellation
%% Simulation settings

timerval = tic;
Delta1 = 0.00001;% precision threshold to solve optimal K_i. unit: second
Delta2 = 10; %search  interval around the value of VectorK_i_T3BH
RepeatTimes = 2000;
RepeatTimes2 = 1000; % the repeat times to find optimal K_i in simulation
Vector_r = 2:4;%the maximum of r should be determined carefully later

alpha0 = 0.5;% unit: second
mu0 = 1/alpha0/10;
N = 200;
L = 20;
d = 4;
KPrimeMax = floor((N - d*(L-1)-1)/d)-1; % to avoid Endless loop  minus  extra 1
VectorKPrime = 2:KPrimeMax; %K should be larger than 2 to make the encoding meaningful'

t_r_Rounding = zeros(length(Vector_r),length(VectorKPrime));
t_r_LB = t_r_Rounding;
t_r_T3BH = t_r_Rounding;
t_r_ED = t_r_Rounding;

DelayLCC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayUncoded = ones(1,length(VectorKPrime)); % the average Delay of the uncoded scheme, which does not change with task segmentation
% DelayAPCC_Optimal = zeros(1,length(VectorKPrime)); % Optimal solution of K, r, K_i, If Without constraints on computation load of each worker

DelayAPCC_r_T3BH = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_SimuOpt = zeros(length(Vector_r),length(VectorKPrime)); % the simulation delay with optimal division (brute force search)

% % If r is too large, the corresponding KPrime, N should also be very large
% % to make (40b) satisfied

DelayAPCCMinimum_T3BH = alpha0 + 1/mu0; % initialize the minimum of delay of APCC for different r
DelayLCCMinimum = alpha0 + 1/mu0; % the minimum of delay of LCC

rOptimalIndex_T3BH = 0;

ExcuTimeLB = 0;
ExcuTimeT3BH = 0;
ExcuTimeSimuOpt = 0;

% DivisionRounding = cell(1,length(Vector_r));
% DivisionLB = cell(1,length(Vector_r));
DivisionT3BH = cell(1,length(Vector_r));
DivisionSimuOpt = cell(1,length(Vector_r));

%% Simulation start

for index_r = 1:length(Vector_r) %the maximum of r should be determined carefully later
    r = Vector_r(index_r);
%     DivisionTest = zeros(r,length(VectorKPrime)); % The division result based on T2 with fixed r
%     DivisionT3BH(1,index_r) = {DivisionTest};
    for index1 = 1:length(VectorKPrime)
        disp(['Progress:',num2str(100*index1/length(VectorKPrime)),'%']);
        KPrime = VectorKPrime(index1);
        EtaPrime = d*(KPrime + L - 1) + 1;
        K = KPrime * r; % the number of division in APCC
        H = d*(K + r*L - r) + r; % according to the relationship between eta_i and K_i in LCC and Case 1 of APCC
        mu = mu0 * K;
        alpha = alpha0 / K; % show the paremeter settings of subtask in APCC

        %% to obtain K_i according to equal division

        VectorK_i_ED = ones(1,r)*K/r;
        VectorK_i_ED = floor(VectorK_i_ED);
        DeltaK = K - sum(VectorK_i_ED); % To ensure the sum of K_i equal to K, basic idea: keep K_i >= K_{i+1}
        while DeltaK ~= 0
            for index2 = 1:r
                if DeltaK == 0
                    break;
                end
                if VectorK_i_ED(index2) < KPrimeMax + 1
                    VectorK_i_ED(index2) = VectorK_i_ED(index2) + 1;
                    DeltaK = DeltaK - 1;
                end
            end
        end
        VectorEta_i_ED = d*(VectorK_i_ED + L -1) +1;

        Vector_t_i = zeros(1,r);% t_i = Vector_t_i(i+1),i is taken within [0:r-1]
        for index3 = 1:r
            Vector_t_i(index3) = alpha*index3 + index3/mu/(N - VectorEta_i_ED(index3) + 1);

        end
        tED = max(Vector_t_i);
        t_r_ED(index_r,index1) = tED;

        %% to obtain rounded K_i according to Theorem 3,
        % which serves as an upperbound for B&B algorithm.

        [VectorKStar_i, Flag, tStar] = T3BasedK_i(H,N,d,L,mu,alpha,r,Delta1);
        % tStar is the lower bound (cannot reach) of B&B algorithm
        % VectorKStar_i can be non-integer, K_i = VectorKStar_i(i+1); Flag indicates whether T2 is solvable
        if Flag == 0
            disp('Unsolvable!')
            continue;
        end
        if abs(K - sum(VectorKStar_i)) > 0.01*K % to make the sum equals K
            disp('Error1!');
        end
        VectorK_i_Rounding = floor(VectorKStar_i); % rounding to obtain integer K_i

        DeltaK = K - sum(VectorK_i_Rounding); % To ensure the sum of K_i equal to K, basic idea: keep K_i >= K_{i+1}
        while DeltaK ~= 0
            for index2 = 1:r
                if DeltaK <= 0
                    break;
                end
                if VectorK_i_Rounding(index2) < KPrimeMax + 1
                    VectorK_i_Rounding(index2) = VectorK_i_Rounding(index2) + 1;
                    DeltaK = DeltaK - 1;
                end
            end
            for index2 = 1:r
                if DeltaK >= 0
                    break;
                end
                if VectorK_i_Rounding(r+1-index2) > 2
                    VectorK_i_Rounding(r+1-index2) = VectorK_i_Rounding(r+1-index2) - 1;
                    DeltaK = DeltaK + 1;
                end
            end
        end
        VectorEta_i_Rounding = d*(VectorK_i_Rounding + L -1) +1; % to update the corresponding integer eta_i, Eta_i = VectorEta_i_Rounding(i+1);
%         DivisionRounding{1,index_r}(:,index1) = VectorK_i_Rounding';

        Vector_t_i = zeros(1,r);% t_i = Vector_t_i(i+1),i is taken within [0:r-1]
        for index3 = 1:r
            Vector_t_i(index3) = alpha*index3 + index3/mu/(N - VectorEta_i_Rounding(index3) + 1);
        end
        tUB = max(Vector_t_i); % serves as an upperbound for B&B algorithm

        t_r_Rounding(index_r,index1) = tUB;

        %% To obtain the optimal integer K_i (LB) using B&B searching algorithm,
        % and the value of corresponding optimization objective serves as the real
        % reachable lowerbound.

        tic
        [VectorK_i_LB, FlagLB, tLB] = T3LowerBoundK_i(K,N,d,L,mu,alpha,r,tUB,r,2,Delta1);
        ExcuTimeLB = ExcuTimeLB + toc;
        VectorEta_i_LB = d*(VectorK_i_LB + L -1) +1;
        t_r_LB(index_r,index1) = tLB;
%         DivisionLB{1,index_r}(:,index1) = VectorK_i_LB';

        %% To obtain the integer K_i (T3BH) using Theorem3 based Heuristic algorithm

        tic
        [VectorK_i_T3BH, FlagT3BH, tT3BH] = T3BH_K_i(K,N,d,L,mu,alpha,r,tUB,r,2,Delta1,VectorKStar_i);
        ExcuTimeT3BH = ExcuTimeT3BH + toc;
        VectorEta_i_T3BH = d*(VectorK_i_T3BH + L -1) +1;
        t_r_T3BH(index_r,index1) = tT3BH;
        DivisionT3BH{1,index_r}(:,index1) = VectorK_i_T3BH';

        %% To obtain the optimal K_i in simulation:APCC with cancellation

        %initilization
        minK_i = max(2,VectorK_i_T3BH(end)-Delta2);
        maxK_i = min([floor(K/r), KPrimeMax + 1, VectorK_i_T3BH(end)+Delta2]);
        VectorK_i = zeros(1,r);
        DelayAPCCMinTemp = alpha0 + 1/mu0;%max(sum(VectorTm));
        VectorK_iMinTemp = VectorK_i_T3BH;
        
%         if r > Vector_r(end)
        if r == Vector_r(end)
            tic
            [DelayAPCCMinTemp, VectorK_iMinTemp] = LoopSearch_w_Cancellation(minK_i,maxK_i,VectorK_i,r,K,DelayAPCCMinTemp,VectorK_iMinTemp,RepeatTimes2,N,d,L,K,mu0,alpha0,VectorK_i_T3BH,Delta2);
            ExcuTimeSimuOpt = ExcuTimeSimuOpt + toc;
        end

        VectorK_i_SimuOpt = VectorK_iMinTemp;
        VectorEta_i_SimuOpt = d*(VectorK_i_SimuOpt + L -1) +1;
        DivisionSimuOpt{1,index_r}(:,index1) = VectorK_i_SimuOpt';


        %% To obtain the expectation of computational delay with cancellation

        SumDelayAPCC_Cancellation = 0;
        SumDelayAPCC_Cancellation_SimuOpt = 0;
        for index3 = 1: RepeatTimes
            u1 = rand(r,N);
            VectorT = alpha0 / K - log(1-u1) / (mu0 * K);
            VectorTm = VectorT;% the delay of each subtask for each worker
            
            DelayAPCCTemp = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i_T3BH); % to obtain the delay of the complete task
            SumDelayAPCC_Cancellation = SumDelayAPCC_Cancellation + DelayAPCCTemp;
            if sum(abs(VectorEta_i_SimuOpt - VectorEta_i_T3BH)) ~= 0 && r == Vector_r(end)
                DelayAPCCTemp2 = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i_SimuOpt);
                SumDelayAPCC_Cancellation_SimuOpt = SumDelayAPCC_Cancellation_SimuOpt + DelayAPCCTemp2;
            end
        end
        DelayAPCC_Cancellation = SumDelayAPCC_Cancellation / RepeatTimes;
        DelayAPCC_Cancellation_SimuOpt = 0;
        if r == Vector_r(end)
            if sum(abs(VectorEta_i_SimuOpt - VectorEta_i_T3BH)) ~= 0
                DelayAPCC_Cancellation_SimuOpt = SumDelayAPCC_Cancellation_SimuOpt / RepeatTimes;
            else
                DelayAPCC_Cancellation_SimuOpt = DelayAPCC_Cancellation;
            end
        end

        DelayAPCC_r_T3BH(index_r,index1) = DelayAPCC_Cancellation;
        DelayAPCC_r_SimuOpt(index_r,index1) = DelayAPCC_Cancellation_SimuOpt;

        if DelayAPCC_Cancellation <= DelayAPCCMinimum_T3BH
            DelayAPCCMinimum_T3BH = DelayAPCC_Cancellation;
            rOptimalIndex_T3BH = index_r;
        end

    end % End of loop for KPrime

    % To obtain the optimal choice of r. Generally speaking, the larger r is,
    % the better.

end% End of loop for r
disp(['Execution time of BF Search:',num2str(ExcuTimeSimuOpt)]);
disp(['Execution time of LB:',num2str(ExcuTimeLB)]);
disp(['Execution time of T3BH:',num2str(ExcuTimeT3BH)]);
DelayAPCC_r_Optimal_T3BH = DelayAPCC_r_T3BH(rOptimalIndex_T3BH,:);

%% Delay of LCC
for index1 = 1:length(VectorKPrime)
    KPrime = VectorKPrime(index1);
    EtaPrime = d*(KPrime + L - 1) + 1;
    SumDelayLCC = 0;
    for index3 = 1: RepeatTimes
        u2 = rand(1,N);
        VectorTPrime = alpha0 / KPrime - log(1-u2) / (mu0 * KPrime);
        [SortVectorTPrime,~] = sort(VectorTPrime);
        DelayLCCTemp = SortVectorTPrime(EtaPrime);
        SumDelayLCC = SumDelayLCC + DelayLCCTemp;
    end
    DelayLCC(index1) = SumDelayLCC / RepeatTimes;
    if DelayLCC(index1) <= DelayLCCMinimum
        DelayLCCMinimum = DelayLCC(index1);
    end
end

%% Delay of Uncoded

% SumDelayUncoded = 0;
% for index3 = 1: RepeatTimes
%     u3 = rand(1,N);
%     VectorTPrime = alpha0 / N - log(1-u3) / (mu0 * N);
%     DelayUncodedTemp = max(VectorTPrime);
%     SumDelayUncoded = SumDelayUncoded + DelayUncodedTemp;
% end
% DelayUncoded = SumDelayUncoded / RepeatTimes *DelayUncoded;

%% Setting of Display figures

Disp_r = Vector_r(rOptimalIndex_T3BH);
t_r_Rounding_Disp = t_r_Rounding(Vector_r == Disp_r,:);
t_r_LB_Disp = t_r_LB(Vector_r == Disp_r,:);
t_r_T3BH_Disp = t_r_T3BH(Vector_r == Disp_r,:);
t_r_ED_Disp = t_r_ED(Vector_r == Disp_r,:);

DelayAPCC1 = DelayAPCC_r_T3BH(Vector_r == 2,:);
DelayAPCC2 = DelayAPCC_r_T3BH(Vector_r == 3,:);
DelayAPCC_SimuOpt1 = DelayAPCC_r_SimuOpt(rOptimalIndex_T3BH,:);

%% Theoretical results

figure(1),hold on

ytemp = 1000*t_r_ED_Disp(t_r_ED_Disp ~= 0);
xtemp = VectorKPrime(t_r_ED_Disp ~= 0);
plot(xtemp,ytemp,'DisplayName',['Equal Division, r=' num2str(Disp_r)])

ytemp = 1000*t_r_Rounding_Disp(t_r_Rounding_Disp ~= 0);
xtemp = VectorKPrime(t_r_Rounding_Disp ~= 0);
plot(xtemp,ytemp,'DisplayName',['T3 Based Rounding, r=' num2str(Disp_r)])

ytemp = 1000*t_r_T3BH_Disp(t_r_T3BH_Disp ~= 0);
xtemp = VectorKPrime(t_r_T3BH_Disp ~= 0);
plot(xtemp,ytemp,'-s','DisplayName',['T3 Based Heuristic, r=' num2str(Disp_r)])

ytemp = 1000*t_r_LB_Disp(t_r_LB_Disp ~= 0);
xtemp = VectorKPrime(t_r_LB_Disp ~= 0);
plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['Lower Bound, r=' num2str(Disp_r)])

line([KPrimeMax+1 KPrimeMax+1],[0 max(t_r_T3BH_Disp)*1200],'linestyle','--','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

legend
title('Theoretical Division Comparison of APCC w/ cancellation')
xlabel('K/r, the number of task segmentation')
ylabel('Value of Objective function: z (ms)')
xlim([-inf KPrimeMax+2])
ylim([0 100])
hold off

%% Simulation results

figure(2),hold on
ytemp = 1000*DelayAPCC1(DelayAPCC1 ~= 0);
xtemp = VectorKPrime(DelayAPCC1 ~= 0);
plot(xtemp,ytemp,'DisplayName',['APCC w/ cancellation, T3BH, r=' num2str(2)])

ytemp = 1000*DelayAPCC2(DelayAPCC2 ~= 0);
xtemp = VectorKPrime(DelayAPCC2 ~= 0);
plot(xtemp,ytemp,'DisplayName',['APCC w/ cancellation, T3BH, r=' num2str(3)])

ytemp = 1000*DelayAPCC_r_Optimal_T3BH(DelayAPCC_r_Optimal_T3BH ~= 0);
xtemp = VectorKPrime(DelayAPCC_r_Optimal_T3BH ~= 0);
plot(xtemp,ytemp,'-s','DisplayName',['APCC w/ cancellation, T3BH, r=' num2str(Vector_r(rOptimalIndex_T3BH))])

ytemp = 1000*DelayAPCC_SimuOpt1(DelayAPCC_SimuOpt1 ~= 0);
xtemp = VectorKPrime(DelayAPCC_SimuOpt1 ~= 0);
plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['APCC w/ cancellation, Optimal(Brute Force), r=' num2str(Vector_r(rOptimalIndex_T3BH))])

% ytemp = DelayUncoded;
% xtemp = VectorKPrime;
% plot(xtemp,ytemp,'LineWidth',1.5,'DisplayName','Uncoded')

ytemp = 1000*DelayLCC;
xtemp = VectorKPrime;
plot(xtemp,ytemp,'-s','DisplayName','LCC')

line([KPrimeMax+1 KPrimeMax+1],[0 max(DelayLCC)*1100],'linestyle','--','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

legend
title('APCC w/ cancellation vs LCC')
xlabel('K''=K/r, the number of task segmentation')
ylabel('Average Delay (ms)')
xlim([-inf KPrimeMax+2])
ylim([0 700])

toc(timerval)

% DelayAPCC_SimuTest_T3BH = zeros(length(Vector_r),length(VectorKPrime));
% Probability_T3BH = zeros(length(Vector_r),length(VectorKPrime));
% DelayAPCC_SimuTest_MaxTimes = zeros(length(Vector_r),length(VectorKPrime));
% Probability_MaxTimes = zeros(length(Vector_r),length(VectorKPrime));

%         DelayAPCCMinExample = alpha0 + 1/mu0;
%         VectorK_iMinExample = zeros(1,r);
%         TimesExample = 1; % the number of times of reaching minimum for some division
%         DelayData = {DelayAPCCMinExample,VectorK_iMinExample,TimesExample};
% 
%         for index3 = 1: RepeatTimes2
%             u1 = rand(r,N);
%             VectorT = alpha0 / K - log(1-u1) / (mu0 * K);
%             VectorTm = VectorT;% the delay of each subtask for each worker

%             DelayDataNum = size(DelayData,1);
%             RepeatMark = 0; % indicates whether there are identical division results
%             for index4 = 1:DelayDataNum
%                 if DelayData{index4,2} == VectorK_iMinTemp
%                     DelayData{index4,1} = (DelayData{index4,1}*DelayData{index4,3} + DelayAPCCMinTemp)/(DelayData{index4,3} + 1);
%                     DelayData{index4,3} = DelayData{index4,3} + 1;
%                     RepeatMark = 1;
%                     break;
%                 end
%             end
%             if RepeatMark == 0
%                 DelayData(DelayDataNum+1,:) = {DelayAPCCMinTemp,VectorK_iMinTemp,1};
%             end
%         end
        
%         DelayDataNum = size(DelayData,1);
%         MaxTimes_index = 1; % the index which corresponds to the division that has max times of reaching delay minimum
%         T3BH_index = 1; % the index corresponding to T3BH division
%         for index5 = 2:DelayDataNum
%             if DelayData{index5,3} >= DelayData{MaxTimes_index,3}
%                 if DelayData{index5,3} > DelayData{MaxTimes_index,3}
%                     MaxTimes_index = index5;
%                 elseif DelayData{index5,1} < DelayData{MaxTimes_index,1}
%                     MaxTimes_index = index5;
%                 end
%             end
%             if DelayData{index5,2} == VectorK_i_T3BH
%                 T3BH_index = index5;
%             end
%         end
%         
%         if T3BH_index ~= 2
%             DelayAPCC_SimuTest_T3BH(index_r,index1) = DelayData{T3BH_index,1};
%             Probability_T3BH(index_r,index1) = DelayData{T3BH_index,3} / RepeatTimes2;
%         else
%             Probability_T3BH(index_r,index1) = 0;
%         end
%         DelayAPCC_SimuTest_MaxTimes(index_r,index1) = DelayData{MaxTimes_index,1};
%         Probability_MaxTimes(index_r,index1) = DelayData{MaxTimes_index,3} / RepeatTimes2;
% 
%         VectorK_i_SimuOpt = DelayData{MaxTimes_index,2};
%         VectorEta_i_SimuOpt = d*(VectorK_i_SimuOpt + L -1) +1;
%         DivisionSimuOpt{1,index_r}(:,index1) = VectorK_i_SimuOpt';