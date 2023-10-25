clc;clear;close all;
%% For comparison with BACC
%% Simulation settings

timerval = tic;
Delta1 = 0.00001;% precision threshold to solve optimal K_i. unit: second
% Delta2 = 5; %search  interval around the value of VectorK_i_T2BH
RepeatTimes = 1000;
% RepeatTimes2 = 1000; % the repeat times to find optimal K_i in simulation
Vector_r = 12;%the maximum of r should be determined carefully later
% Realtionship between eta_i and K_i: L=25 d = 6 for d =25，L=0, N=500, with error
%controlled in 5%, 
%L = 35 d = 6 for N =700 d =35

alpha0 = 0.5;% unit: second
mu0 = 1/alpha0/10;
N = 500;
L = 50;
d = 4;% d =25
% KPrimeMax = floor(N/12)-1;
KPrimeMax = floor((N - d*(L-1)-1)/d)-1; % to avoid Endless loop  minus  extra 1
VectorKPrime = 2:KPrimeMax; %K should be larger than 2 to make the encoding meaningful'

t_r_Rounding = zeros(length(Vector_r),length(VectorKPrime));
t_r_T2BH = t_r_Rounding;
t_r_T3BH = t_r_Rounding;

DelayLCC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayBACC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayUncoded = ones(1,length(VectorKPrime)); % the average Delay of the uncoded scheme, which does not change with task segmentation
% DelayAPCC_Optimal = zeros(1,length(VectorKPrime)); % Optimal solution of K, r, K_i, If Without constraints on computation load of each worker

DelayAPCC_r_T2BH = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_T3BH = zeros(length(Vector_r),length(VectorKPrime));
% DelayAPCC_r_SimuOpt = zeros(length(Vector_r),length(VectorKPrime)); % the simulation delay with optimal division (brute force search)

% % If r is too large, the corresponding KPrime, N should also be very large
% % to make (40b) satisfied

DelayAPCCMinimum_T2BH = alpha0 + 1/mu0; % initialize the minimum of delay of APCC for different r
DelayAPCCMinimum_T3BH = alpha0 + 1/mu0; % initialize the minimum of delay of APCC for different r
DelayLCCMinimum = alpha0 + 1/mu0; % the minimum of delay of LCC
DelayBACCMinimum = alpha0 + 1/mu0; % the minimum of delay of BACC

rOptimalIndex_T2BH = 0;
rOptimalIndex_T3BH = 0;

ExcuTimeT2BH = 0;
ExcuTimeT3BH = 0;

DivisionT2BH = cell(1,length(Vector_r));
DivisionT3BH = cell(1,length(Vector_r));

%% Simulation start

for index_r = 1:length(Vector_r) %the maximum of r should be determined carefully later
    r = Vector_r(index_r);
    %     DivisionTest = zeros(r,length(VectorKPrime)); % The division result based on T2 with fixed r
    for index1 = 1:length(VectorKPrime)
        disp(['Progress:',num2str(100*index1/length(VectorKPrime)),'%']);
        KPrime = VectorKPrime(index1);
        EtaPrime = d*(KPrime + L - 1) + 1;
%         EtaPrime = 12*KPrime;
        K = KPrime * r; % the number of division in APCC
%         H = 12*K;
        H = d*(K + r*L - r) + r; % according to the relationship between eta_i and K_i in LCC and Case 1 of APCC
        mu = mu0 * K;
        alpha = alpha0 / K; % show the paremeter settings of subtask in APCC

        %% to obtain rounded K_i according to Theorem 2,
        % which serves as an upperbound for B&B algorithm.
        
%         if r==5 && K == 15
%             disp('test');
%         end
            
        [VectorKStar_i, Flag, ~] = T2BasedK_i(H,N,d,L,mu,alpha,r,Delta1);
        % tStar is the lower bound (cannot reach) of B&B algorithm
        % VectorKStar_i can be non-integer, K_i = VectorKStar_i(i+1); Flag indicates whether T2 is solvable
        if Flag == 0
            disp('Unsolvable!')
            VectorKStar_i = ones(1,r)*K/r; %if unsolvable, replace with the equal division scheme
        end
        if abs(K - sum(VectorKStar_i)) > 0.01*K % to make the sum equals K
            disp('Error1!');
        end
        VectorK_i_Rounding = floor(VectorKStar_i); % rounding to obtain integer K_i

        VectorK_i_Rounding = RoundingK_i(K,r,VectorK_i_Rounding,KPrimeMax);
        VectorEta_i_Rounding = d*(VectorK_i_Rounding + L -1) +1; % to update the corresponding integer eta_i, Eta_i = VectorEta_i_Rounding(i+1);
        %         DivisionRounding{1,index_r}(:,index1) = VectorK_i_Rounding';

        Vector_t_i = zeros(1,r);% t_i = Vector_t_i(i+1),i is taken within [0:r-1]
        for index3 = 1:r
            Vector_t_i(index3) = alpha*index3 - index3*log(1 - VectorEta_i_Rounding(index3)/N)/mu;
        end
        tUB = max(Vector_t_i); % serves as an upperbound for B&B algorithm

        t_r_Rounding(index_r,index1) = tUB;

        %% To obtain the integer K_i (T2BH) using Theorem2 based Heuristic algorithm

        tic
        [VectorK_i_T2BH, FlagT2BH, tT2BH] = T2BH_K_i(K,N,d,L,mu,alpha,r,tUB,VectorK_i_Rounding,r,2,Delta1,VectorKStar_i);
        %         disp(['Execution time of LB:',num2str(toc)]);
        ExcuTimeT2BH = ExcuTimeT2BH + toc;
        VectorEta_i_T2BH = d*(VectorK_i_T2BH + L -1) +1;
        DivisionT2BH{1,index_r}(:,index1) = VectorK_i_T2BH';
        t_r_T2BH(index_r,index1) = tT2BH;

        %% to obtain rounded K_i according to Theorem 3,
        % which serves as an upperbound for B&B algorithm.

        [VectorKStar_i, Flag, tStar] = T3BasedK_i(H,N,d,L,mu,alpha,r,Delta1);
        % tStar is the lower bound (cannot reach) of B&B algorithm
        % VectorKStar_i can be non-integer, K_i = VectorKStar_i(i+1); Flag indicates whether T2 is solvable
        if Flag == 0
            disp('Unsolvable!')
            VectorKStar_i = ones(1,r)*K/r; %if unsolvable, replace with the equal division scheme
        end
        if abs(K - sum(VectorKStar_i)) > 0.01*K % to make the sum equals K
            disp('Error1!');
        end
        VectorK_i_Rounding = floor(VectorKStar_i); % rounding to obtain integer K_i

        VectorK_i_Rounding = RoundingK_i(K,r,VectorK_i_Rounding,KPrimeMax);
        VectorEta_i_Rounding = d*(VectorK_i_Rounding + L -1) +1; % to update the corresponding integer eta_i, Eta_i = VectorEta_i_Rounding(i+1);
        %         DivisionRounding{1,index_r}(:,index1) = VectorK_i_Rounding';

        Vector_t_i = zeros(1,r);% t_i = Vector_t_i(i+1),i is taken within [0:r-1]
        for index3 = 1:r
            Vector_t_i(index3) = alpha*index3 + index3/mu/(N - VectorEta_i_Rounding(index3) + 1);
        end
        tUB = max(Vector_t_i); % serves as an upperbound for B&B algorithm

        t_r_Rounding(index_r,index1) = tUB;

        %% To obtain the integer K_i (T3BH) using Theorem3 based Heuristic algorithm

        tic
        [VectorK_i_T3BH, FlagT3BH, tT3BH] = T3BH_K_i(K,N,d,L,mu,alpha,r,tUB,VectorK_i_Rounding,r,2,Delta1,VectorKStar_i);
        ExcuTimeT3BH = ExcuTimeT3BH + toc;
        VectorEta_i_T3BH = d*(VectorK_i_T3BH + L -1) +1;
        t_r_T3BH(index_r,index1) = tT3BH;
        DivisionT3BH{1,index_r}(:,index1) = VectorK_i_T3BH';

        %% To obtain the expectation of computational delay with/without cancellation

        SumDelayAPCC = 0;
        SumDelayAPCC_Cancellation = 0;
        for index3 = 1: RepeatTimes
            u1 = rand(r,N);
            VectorT = alpha0 / K - log(1-u1) / (mu0 * K);
            VectorTm = VectorT;

            DelayAPCCTemp = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i_T2BH); % to obtain the delay of the complete task
            SumDelayAPCC = SumDelayAPCC + DelayAPCCTemp;
            DelayAPCCTemp2 = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i_T3BH);
            SumDelayAPCC_Cancellation = SumDelayAPCC_Cancellation + DelayAPCCTemp2;

        end
        DelayAPCC = SumDelayAPCC / RepeatTimes;
        DelayAPCC_Cancellation = SumDelayAPCC_Cancellation / RepeatTimes;

        DelayAPCC_r_T2BH(index_r,index1) = DelayAPCC;
        DelayAPCC_r_T3BH(index_r,index1) = DelayAPCC_Cancellation;

        if DelayAPCC <= DelayAPCCMinimum_T2BH
            DelayAPCCMinimum_T2BH = DelayAPCC;
            rOptimalIndex_T2BH = index_r;
        end
        if DelayAPCC_Cancellation <= DelayAPCCMinimum_T3BH
            DelayAPCCMinimum_T3BH = DelayAPCC_Cancellation;
            rOptimalIndex_T3BH = index_r;
        end

    end % End of loop for KPrime
    % To obtain the optimal choice of r. Generally speaking, the larger r is,
    % the better.
end% End of loop for r
disp(['Execution time of T2BH:',num2str(ExcuTimeT2BH)]);
disp(['Execution time of T3BH:',num2str(ExcuTimeT3BH)]);
DelayAPCC_r_Optimal_T2BH = DelayAPCC_r_T2BH(rOptimalIndex_T2BH,:);
DelayAPCC_r_Optimal_T3BH = DelayAPCC_r_T3BH(rOptimalIndex_T3BH,:);

%% Delay of Uncoded

SumDelayUncoded = 0;
for index3 = 1: RepeatTimes
    u3 = rand(1,N);
    VectorTPrime = alpha0 / N - log(1-u3) / (mu0 * N);
    DelayUncodedTemp = max(VectorTPrime);
    SumDelayUncoded = SumDelayUncoded + DelayUncodedTemp;
end
DelayUncoded = SumDelayUncoded / RepeatTimes *DelayUncoded;

%% Delay of BACC

for index1 = 1:length(VectorKPrime)
    KPrime = VectorKPrime(index1);
    %     EtaPrime = d/4*(KPrime + L );
    %     EtaPrime = 4*(KPrime + L );
    EtaPrime = d*(KPrime + L - 1) + 1;
    SumDelayBACC = 0;
    for index3 = 1: RepeatTimes
        u2 = rand(1,N);
        VectorTPrime = alpha0 / KPrime - log(1-u2) / (mu0 * KPrime);
        [SortVectorTPrime,~] = sort(VectorTPrime);
        DelayBACCTemp = SortVectorTPrime(EtaPrime);
        SumDelayBACC = SumDelayBACC + DelayBACCTemp;
    end
    DelayBACC(index1) = SumDelayBACC / RepeatTimes;
    if DelayBACC(index1) <= DelayBACCMinimum
        DelayBACCMinimum = DelayBACC(index1);
    end
end

%% Display figures

DelayAPCC1 = DelayAPCC_r_T2BH(Vector_r == 2,:);
DelayAPCC2 = DelayAPCC_r_T3BH(Vector_r == 2,:);


%% Simulation results

figure(1),hold on
% ytemp = 1000*DelayAPCC1(DelayAPCC1 ~= 0);
% xtemp = VectorKPrime(DelayAPCC1 ~= 0);
% plot(xtemp,ytemp,'DisplayName',['APCC w/o cancellation, T2BH, r=' num2str(2)])
% 
% ytemp = 1000*DelayAPCC2(DelayAPCC2 ~= 0);
% xtemp = VectorKPrime(DelayAPCC2 ~= 0);
% plot(xtemp,ytemp,'DisplayName',['APCC w/ cancellation, T2BH, r=' num2str(2)])

ytemp = 1000*DelayAPCC_r_Optimal_T2BH(DelayAPCC_r_Optimal_T2BH ~= 0);
xtemp = VectorKPrime(DelayAPCC_r_Optimal_T2BH ~= 0);
plot(xtemp,ytemp,'DisplayName',['APCC w/o cancellation, T2BH, r=' num2str(Vector_r(rOptimalIndex_T2BH))])

ytemp = 1000*DelayAPCC_r_Optimal_T3BH(DelayAPCC_r_Optimal_T3BH ~= 0);
xtemp = VectorKPrime(DelayAPCC_r_Optimal_T3BH ~= 0);
plot(xtemp,ytemp,'-','DisplayName',['APCC w/ cancellation, T3BH, r=' num2str(Vector_r(rOptimalIndex_T3BH))])

ytemp = 1000*DelayUncoded;
xtemp = VectorKPrime;
plot(xtemp,ytemp,'DisplayName','Uncoded')

ytemp = 1000*DelayBACC;
xtemp = VectorKPrime;
plot(xtemp,ytemp,'-','DisplayName','BACC')

line([KPrimeMax+1 KPrimeMax+1],[0 max(DelayBACC)*1100],'linestyle','-','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

legend
title('APCC vs BACC')
xlabel('K''=K/r, the number of task segmentation')
ylabel('Average Delay (ms)')
xlim([-inf KPrimeMax+2])
% ylim([0 max(DelayBACC)*1100])
ylim([0 400])

toc(timerval)