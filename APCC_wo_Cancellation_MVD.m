clc;clear;close all;
%% For performance verification of APCC without cancellation
%% Simulation settings

timerval = tic;
Delta1 = 0.00001;% precision threshold to solve optimal K_i. unit: second
Delta2 = 5; %search  interval around the value of VectorK_i_T2BH
Mark_TheoryResult = 0;
Mark_simulation = 1; %decide whether to conduct simulation, 1:yes 0:no
Mark_LB = 0; %decide whether to perform LB algorithm.
Mark_T2BH = 0;
Mark_BF = 0; 
RepeatTimes = 2000;
RepeatTimes2 = 1000; % the repeat times to find optimal K_i in simulation
Vector_r = [4,8,16];%the maximum of r should be determined carefully later
r_BF = 4; % the r that corresponds to APCC,BF algorithm

alpha0 = 0.5;% unit: second
mu0 = 1/alpha0/10;
N = 200;
L = 20;
d = 4;
KPrimeMax = floor((N - d*(L-1)-1)/d)-1; % to avoid Endless loop  minus  extra 1
VectorKPrime = 2:KPrimeMax; %K should be larger than 2 to make the encoding meaningful'

t_r_Rounding = zeros(length(Vector_r),length(VectorKPrime));
t_r_LB = t_r_Rounding;
t_r_T2BH = t_r_Rounding;
t_r_MVD = t_r_Rounding;
t_r_ED = t_r_Rounding;

DelayLCC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayUncoded = ones(1,length(VectorKPrime)); % the average Delay of the uncoded scheme, which does not change with task segmentation
% DelayAPCC_Optimal = zeros(1,length(VectorKPrime)); % Optimal solution of K, r, K_i, If Without constraints on computation load of each worker

DelayAPCC_r_MVD = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_SimuOpt = zeros(length(Vector_r),length(VectorKPrime)); % the simulation delay with optimal division (brute force search)

% % If r is too large, the corresponding KPrime, N should also be very large
% % to make (40b) satisfied

 % initialize the minimum of delay of APCC for different r
DelayLCCMinimum = alpha0 + 1/mu0; % the minimum of delay of LCC
% DelayAPCCMinimum_MVD = alpha0 + 1/mu0;
Vec_r_DelayAPCCMinimum_MVD = ones(1,length(Vector_r))*(alpha0 + 1/mu0);
Vec_r_DelayAPCCMinimum_SimuOpt = ones(1,length(Vector_r))*(alpha0 + 1/mu0);

% r_Optimal_MVD = 0;

ExcuTimeLB = 0;
ExcuTimeT2BH = 0;
ExcuTimeSimuOpt = 0;
ExcuTimeMVD = 0;

% DivisionRounding = cell(1,length(Vector_r));
DivisionLB = cell(1,length(Vector_r));
% DivisionT2BH = cell(1,length(Vector_r));
DivisionMVD = cell(1,length(Vector_r));
DivisionSimuOpt = cell(1,length(Vector_r));

%% Simulation start

for index_r = 1:length(Vector_r) %the maximum of r should be determined carefully later
    r = Vector_r(index_r);
    VectorK = r*VectorKPrime;
    %     DivisionTest = zeros(r,length(VectorKPrime)); % The division result based on T2 with fixed r
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
            Vector_t_i(index3) = alpha*index3 - index3*log(1 - VectorEta_i_ED(index3)/N)/mu;
        end
        tED = max(Vector_t_i);
        t_r_ED(index_r,index1) = tED;

        %% to obtain rounded K_i according to Theorem 2,
        % which serves as an upperbound for B&B algorithm.

        [VectorKStar_i, Flag, tStar] = T2BasedK_i(H,N,d,L,mu,alpha,r,Delta1);
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
            Vector_t_i(index3) = alpha*index3 - index3*log(1 - VectorEta_i_Rounding(index3)/N)/mu;
        end
        tUB = max(Vector_t_i); % serves as an upperbound for B&B algorithm

        t_r_Rounding(index_r,index1) = tUB;

        %% To obtain the optimal integer K_i (LB) using B&B searching algorithm,
        % and the value of optimization objective serves as the real
        % reachable lowerbound.

        if Mark_LB == 1
            tic
            [VectorK_i_LB, FlagLB, tLB] = T2LowerBoundK_i(K,N,d,L,mu,alpha,r,tUB,r,2,Delta1);
            %         disp(['Execution time of LB:',num2str(toc)]);
            ExcuTimeLB = ExcuTimeLB + toc;
            VectorEta_i_LB = d*(VectorK_i_LB + L -1) +1;
            DivisionLB{1,index_r}(:,index1) = VectorK_i_LB';
            t_r_LB(index_r,index1) = tLB;
        end

        %% To obtain the integer K_i (T2BH) using Theorem2 based Heuristic algorithm

        if Mark_T2BH == 1
            tic
            [VectorK_i_T2BH, FlagT2BH, tT2BH] = T2BH_K_i(K,N,d,L,mu,alpha,r,tUB,VectorK_i_Rounding,r,2,Delta1,VectorKStar_i);
            %         disp(['Execution time of LB:',num2str(toc)]);
            ExcuTimeT2BH = ExcuTimeT2BH + toc;
            VectorEta_i_T2BH = d*(VectorK_i_T2BH + L -1) +1;
            %         DivisionT2BH{1,index_r}(:,index1) = VectorK_i_T2BH';
            t_r_T2BH(index_r,index1) = tT2BH;
        end

        %% To obtain the integer K_i (MVD) using Maximum Value Descent Algorithm

        tic
        [VectorK_i_MVD, FlagMVD, tMVD] = MVD_T2_K_i(VectorK_i_Rounding,r,K,N,d,L,alpha,mu,KPrimeMax);
        %         disp(['Execution time of LB:',num2str(toc)]);
        ExcuTimeMVD = ExcuTimeMVD + toc;
        VectorEta_i_MVD = d*(VectorK_i_MVD + L -1) +1;
        DivisionMVD{1,index_r}(:,index1) = VectorK_i_MVD';
        t_r_MVD(index_r,index1) = tMVD;

%         if tMVD ~= tLB
% %         if sum(abs(VectorK_i_MVD - VectorK_i_LB)) ~= 0 % this is normal because sometimes the same zStar corresponds to different Vector_K_i
%             disp('MVD difference!'); 
%         end

        %% To obtain the optimal K_i (SimuOpt) in simulation:APCC without cancellation

        %initilization
        minK_i = max(2,VectorK_i_MVD(end)-Delta2);
        maxK_i = min([floor(K/r), KPrimeMax + 1, VectorK_i_MVD(end)+Delta2]);
        VectorK_i = zeros(1,r);
        DelayAPCCMinTemp = alpha0 + 1/mu0;%max(sum(VectorTm));
        VectorK_iMinTemp = VectorK_i_MVD;

%       % decide whether to run this part of codes
        if Mark_BF == 1 && r == r_BF
            tic
            [DelayAPCCMinTemp, VectorK_iMinTemp] = LoopSearch_wo_Cancellation(minK_i,maxK_i,VectorK_i,r,K,DelayAPCCMinTemp,VectorK_iMinTemp,RepeatTimes2,N,d,L,K,mu0,alpha0,VectorK_i_MVD,Delta2);
            ExcuTimeSimuOpt = ExcuTimeSimuOpt + toc;
        end

        VectorK_i_SimuOpt = VectorK_iMinTemp;
        VectorEta_i_SimuOpt = d*(VectorK_i_SimuOpt + L -1) +1;
        DivisionSimuOpt{1,index_r}(:,index1) = VectorK_i_SimuOpt';

        %% To obtain the expectation of computational delay without cancellation
        if Mark_simulation == 1

            %         SumDelayAPCC_Rounding = 0;
            %         SumDelayAPCC_LB = 0;
            SumDelayAPCC = 0;
            SumDelayAPCC_SimuOpt = 0;
            for index3 = 1: RepeatTimes
                u1 = rand(r,N);
                VectorT = alpha0 / K - log(1-u1) / (mu0 * K);
                VectorTm = VectorT;

                DelayAPCCTemp = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i_MVD); % to obtain the delay of the complete task
                SumDelayAPCC = SumDelayAPCC + DelayAPCCTemp;
                if sum(abs(VectorEta_i_SimuOpt - VectorEta_i_MVD)) ~= 0 && r == r_BF
                    DelayAPCCTemp2 = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i_SimuOpt);
                    SumDelayAPCC_SimuOpt = SumDelayAPCC_SimuOpt + DelayAPCCTemp2;
                end
            end
            DelayAPCC = SumDelayAPCC / RepeatTimes;
            if sum(abs(VectorEta_i_SimuOpt - VectorEta_i_MVD)) ~= 0 && r == r_BF
                DelayAPCC_SimuOpt = SumDelayAPCC_SimuOpt / RepeatTimes;
            else
                DelayAPCC_SimuOpt = DelayAPCC;
            end

            DelayAPCC_r_MVD(index_r,index1) = DelayAPCC;
            DelayAPCC_r_SimuOpt(index_r,index1) = DelayAPCC_SimuOpt;

            if DelayAPCC <= Vec_r_DelayAPCCMinimum_MVD(index_r)
                Vec_r_DelayAPCCMinimum_MVD(index_r) = DelayAPCC;
            end
            if DelayAPCC_SimuOpt <= Vec_r_DelayAPCCMinimum_SimuOpt(index_r)
                Vec_r_DelayAPCCMinimum_SimuOpt(index_r) = DelayAPCC_SimuOpt;
            end
        end

    end % End of loop for KPrime

    % To obtain the optimal choice of r. Generally speaking, the larger r is,
    % the better.

end% End of loop for r
disp(['Execution time of BF Search:',num2str(ExcuTimeSimuOpt)]);
disp(['Execution time of LB:',num2str(ExcuTimeLB)]);
disp(['Execution time of T2BH:',num2str(ExcuTimeT2BH)]);
disp(['Execution time of MVD:',num2str(ExcuTimeMVD)]);
[DelayAPCCMinimum_MVD,temp2] = min(Vec_r_DelayAPCCMinimum_MVD);
r_Optimal_MVD = Vector_r(temp2);
t_r_MVD_Min = min(t_r_MVD(end,:));

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

%% Setting of Display figures

Disp_r = Vector_r(end);
t_r_Rounding_Disp = t_r_Rounding(Vector_r == Disp_r,:);
t_r_LB_Disp = t_r_LB(Vector_r == Disp_r,:);
t_r_T2BH_Disp = t_r_T2BH(Vector_r == Disp_r,:);
t_r_MVD_Disp = t_r_MVD(Vector_r == Disp_r,:);
t_r_ED_Disp = t_r_ED(Vector_r == Disp_r,:);

%% Theoretical results

if Mark_TheoryResult == 1

    figure(1),hold on

    ytemp = 1000*t_r_ED_Disp(t_r_ED_Disp ~= 0);
    xtemp = VectorK(t_r_ED_Disp ~= 0);
    plot(xtemp,ytemp,'DisplayName',['Equal Division, r=' num2str(Disp_r)])

    ytemp = 1000*t_r_MVD_Disp(t_r_MVD_Disp ~= 0);
    xtemp = VectorK(t_r_MVD_Disp ~= 0);
    plot(xtemp,ytemp,'-s','DisplayName',['MVD, r=' num2str(Disp_r)])

    ytemp = 1000*t_r_LB_Disp(t_r_LB_Disp ~= 0);
    xtemp = VectorK(t_r_LB_Disp ~= 0);
    plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['B&B, r=' num2str(Disp_r)])

    % line([KPrimeMax+1 KPrimeMax+1],[0 max(t_r_ED_Disp)*1100],'linestyle','--','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

    legend
    title('Theoretical Division Comparison of APCC w/o cancellation')
    xlabel('K, the number of task segmentation')
    ylabel('Value of Objective function: z (ms)')
    xlim([-inf r*(KPrimeMax+2)])
    % ylim([0 700])
    hold off

end

%% Simulation results

if Mark_simulation == 1

    Disp_r1 = 4;
    Disp_r2 = 8;
    DelayAPCC1 = DelayAPCC_r_MVD(Vector_r == Disp_r1,:);
    DelayAPCC2 = DelayAPCC_r_MVD(Vector_r == Disp_r2,:);
    DelayAPCC_r_Optimal_MVD = DelayAPCC_r_MVD(Vector_r == r_Optimal_MVD,:);
    DelayAPCC_SimuOpt1 = DelayAPCC_r_SimuOpt(Vector_r == r_BF,:);


    figure(2),hold on
    ytemp = 1000*DelayAPCC1(DelayAPCC1 ~= 0);
    xtemp = VectorKPrime(DelayAPCC1 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/o cancellation, MVD, r=' num2str(Disp_r1)])

    ytemp = 1000*DelayAPCC2(DelayAPCC2 ~= 0);
    xtemp = VectorKPrime(DelayAPCC2 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/o cancellation, MVD, r=' num2str(Disp_r2)])

    ytemp = 1000*DelayAPCC_r_Optimal_MVD(DelayAPCC_r_Optimal_MVD ~= 0);
    xtemp = VectorKPrime(DelayAPCC_r_Optimal_MVD ~= 0);
    plot(xtemp,ytemp,'-','LineWidth',1.25,'DisplayName',['APCC w/o cancellation, MVD, r=' num2str(r_Optimal_MVD)])

%     ytemp = 1000*DelayAPCC_SimuOpt1(DelayAPCC_SimuOpt1 ~= 0);
%     xtemp = VectorKPrime(DelayAPCC_SimuOpt1 ~= 0);
%     plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['APCC w/o cancellation, Optimal(Brute Force), r=' num2str(r_Optimal_MVD)])

    % ytemp = DelayUncoded;
    % xtemp = VectorKPrime;
    % plot(xtemp,ytemp,'LineWidth',1.5,'DisplayName','Uncoded')

    ytemp = 1000*DelayLCC;
    xtemp = VectorKPrime;
    plot(xtemp,ytemp,'-','LineWidth',1.25,'DisplayName','LCC')

%     line([KPrimeMax+1 KPrimeMax+1],[0 max(DelayLCC)*1100],'linestyle','--','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

    legend
%     title('APCC w/o cancellation vs LCC')
    xlabel('K''=K/r, the number of task segmentation')
    ylabel('Average Delay (ms)')
    xlim([0 KPrimeMax+1])
    ylim([0 800])

end

toc(timerval)