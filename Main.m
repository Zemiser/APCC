clc;clear; close all;
%% For comparison with LCC and LCC-MMC updated on 2023-10-16
%% Simulation settings

timerval = tic;
Delta1 = 0.00001;% precision threshold to solve optimal K_i. unit: second
Delta2 = 5;%round(0.5*(KPrimeMax + 1)); %search  interval around the value of VectorK_i_MVD_T2
Delta3 = 8; %search  interval around the value of VectorK_i_MVD_T3
Mark_simulation = 1;
Mark_Histogram = 1;
Mark_BF = 0; %decide whether to perform Brute-force algorithm.
Mark_rL = 0; %decide whether to perform the analysis of r and L.

RepeatTimes = 2000;
RepeatTimes2 = 1000; % the repeat times to find optimal K_i in simulation
Vector_r = [4,16];%the maximum of r should be determined carefully later
r_BF = 4; % the r that corresponds to APCC,BF algorithm
r_disp = Vector_r(end);

alpha0 = 0.5;% unit: second
mu0 = 1/alpha0/10;
N = 200;
L = 0;
d = 4;
KPrimeMax = floor((N - d*(L-1)-1)/d)-1; % to avoid Endless loop  minus  extra 1
VectorKPrime = 2:KPrimeMax+1; %K should be larger than 2 to make the encoding meaningful'

t_r_MVD_T2 = zeros(length(Vector_r),length(VectorKPrime));
t_r_MVD_T3 = t_r_MVD_T2;

DelayLCC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayLCC_MMC = zeros(length(Vector_r),length(VectorKPrime));
DelayBACC = zeros(1,length(VectorKPrime)); % the average Delay, the same for APCC
DelayUncoded = ones(1,length(VectorKPrime)); % the average Delay of the uncoded scheme, which does not change with task segmentation
% DelayAPCC_Optimal = zeros(1,length(VectorKPrime)); % Optimal solution of K, r, K_i, If Without constraints on computation load of each worker

DelayAPCC_r_MVD_T2 = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_MVD_T3 = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_SimuOpt_T2 = zeros(length(Vector_r),length(VectorKPrime));
DelayAPCC_r_SimuOpt_T3 = zeros(length(Vector_r),length(VectorKPrime)); % the simulation delay with optimal division (brute force search)

% % If r is too large, the corresponding KPrime, N should also be very large
% % to make (40b) satisfied

% DelayLCCMinimum = alpha0 + 1/mu0; % the minimum of delay of LCC
DelayBACCMinimum = alpha0 + 1/mu0; % the minimum of delay of BACC
Vec_r_DelayAPCCMinimum_MVD_T2 = ones(1,length(Vector_r))*(alpha0 + 1/mu0);
Vec_r_DelayAPCCMinimum_SimuOpt_T2 = ones(1,length(Vector_r))*(alpha0 + 1/mu0);
Vec_r_DelayAPCCMinimum_MVD_T3 = ones(1,length(Vector_r))*(alpha0 + 1/mu0);
Vec_r_DelayAPCCMinimum_SimuOpt_T3 = ones(1,length(Vector_r))*(alpha0 + 1/mu0);
Vec_r_DelayLCC_MMC_Minimum = ones(1,length(Vector_r))*(alpha0 + 1/mu0);

rOptimalIndex_MVD_T2 = 0;
rOptimalIndex_MVD_T3 = 0;

ExcuTimeMVD_T2 = 0;
ExcuTimeMVD_T3 = 0;
ExcuTimeSimuOpt_T3 = 0;
ExcuTimeSimuOpt_T2 = 0;

DivisionMVD_T2 = cell(1,length(Vector_r));
DivisionMVD_T3 = cell(1,length(Vector_r));
DivisionSimuOpt_T2 = cell(1,length(Vector_r));
DivisionSimuOpt_T3 = cell(1,length(Vector_r));

%% Simulation start

for index_r = 1:length(Vector_r) %the maximum of r should be determined carefully later

    r = Vector_r(index_r);
    %     DivisionTest = zeros(r,length(VectorKPrime)); % The division result based on T2 with fixed r

    for index1 = 1:length(VectorKPrime)

        disp(['Progress:',num2str(100*index1/length(VectorKPrime)),'%']);
        KPrime = VectorKPrime(index1);
        EtaPrime = d*(KPrime + L - 1) + 1;
        K = KPrime * r; % the number of division in APCC
        K_LM = K; % the number of division in LLC-MMC
        H = d*(K + r*L - r) + r; % according to the relationship between eta_i and K_i in LCC and Case 1 of APCC
        mu = mu0 * K;
        alpha = alpha0 / K; % show the paremeter settings of subtask in APCC

        %% to obtain rounded K_i according to Theorem 2,
        % which serves as an upperbound for B&B algorithm.

        [VectorKStar_i, Flag, ~] = T2BasedK_i(H,N,d,L,mu,alpha,r,Delta1);

        % VectorKStar_i can be non-integer, K_i = VectorKStar_i(i+1); Flag indicates whether T2 is solvable
        if Flag == 0
            disp('Unsolvable!')
            pause
            VectorKStar_i = ones(1,r)*K/r; %if unsolvable, replace with the equal division scheme
        end
        if abs(K - sum(VectorKStar_i)) > 0.01*K % to make the sum equals K
            disp('Error1!');
            pause
        end
        VectorK_i_Rounding = round(VectorKStar_i); % rounding to obtain integer K_i

        VectorK_i_Rounding = RoundingK_i(K,r,VectorK_i_Rounding,KPrimeMax);

        %% To obtain the integer K_i (MVD_T2) using Maximum Value Descent Algorithm

        tic
        [VectorK_i_MVD_T2, ~, tMVD_T2] = MVD_T2_K_i(VectorK_i_Rounding,r,K,N,d,L,alpha,mu,KPrimeMax);
        %         disp(['Execution time of LB:',num2str(toc)]);
        ExcuTimeMVD_T2 = ExcuTimeMVD_T2 + toc;
        VectorEta_i_MVD_T2 = d*(VectorK_i_MVD_T2 + L -1) +1;
        DivisionMVD_T2{1,index_r}(:,index1) = VectorK_i_MVD_T2';
        t_r_MVD_T2(index_r,index1) = tMVD_T2;

        %% To obtain the optimal K_i (SimuOpt_T2) in simulation:APCC without cancellation

        %       % decide whether to run this part of codes
        if Mark_BF == 1 && r == r_BF
            tic
            [K_simu_opt_vec, t_simu_opt] = Brute_force_wo_Cancel(VectorK_i_MVD_T2,r,K,N,d,L,...
                mu0,alpha0,Delta2,KPrimeMax,RepeatTimes2);
            ExcuTimeSimuOpt_T2 = ExcuTimeSimuOpt_T2 + toc;
        else
            K_simu_opt_vec = VectorK_i_MVD_T2;
        end

        VectorK_i_SimuOpt_T2 = K_simu_opt_vec;
        VectorEta_i_SimuOpt_T2 = d*(VectorK_i_SimuOpt_T2 + L -1) +1;
        DivisionSimuOpt_T2{1,index_r}(:,index1) = VectorK_i_SimuOpt_T2';


        %% to obtain rounded K_i according to Theorem 3,
        % which serves as an upperbound for B&B algorithm.

        [VectorKStar_i_T3, Flag, ~] = T3BasedK_i(H,N,d,L,mu,alpha,r,Delta1);

        % VectorKStar_i can be non-integer, K_i = VectorKStar_i(i+1); Flag indicates whether T2 is solvable
        if Flag == 0
            disp('Unsolvable!')
            pause
            VectorKStar_i_T3 = ones(1,r)*K/r; %if unsolvable, replace with the equal division scheme
        end
        if abs(K - sum(VectorKStar_i_T3)) > 0.01*K % to make the sum equals K
            disp('Error1!');
            pause
        end
        VectorK_i_Rounding_T3 = round(VectorKStar_i_T3); % rounding to obtain integer K_i

        VectorK_i_Rounding_T3 = RoundingK_i(K,r,VectorK_i_Rounding_T3,KPrimeMax);

        %% To obtain the integer K_i (MVD_T3) using Maximum Value Descent Algorithm

        tic
        [VectorK_i_MVD_T3, ~, tMVD_T3] = MVD_T3_K_i(VectorK_i_Rounding_T3,r,K,N,d,L,alpha,mu,KPrimeMax);
        ExcuTimeMVD_T3 = ExcuTimeMVD_T3 + toc;
        VectorEta_i_MVD_T3 = d*(VectorK_i_MVD_T3 + L -1) +1;
        t_r_MVD_T3(index_r,index1) = tMVD_T3;
        DivisionMVD_T3{1,index_r}(:,index1) = VectorK_i_MVD_T3';

        %% To obtain the optimal K_i (SimuOpt_T3) in simulation:APCC with cancellation

        %       % decide whether to run this part of codes
        if Mark_BF == 1 && r == r_BF
            tic
            [K_simu_opt_vec, t_simu_opt] = Brute_force_w_Cancel(VectorK_i_MVD_T3,r,K,N,d,L,...
                mu0,alpha0,Delta3,KPrimeMax,RepeatTimes2);
            ExcuTimeSimuOpt_T3 = ExcuTimeSimuOpt_T3 + toc;
        else
            K_simu_opt_vec = VectorK_i_MVD_T3;
        end

        VectorK_i_SimuOpt_T3 = K_simu_opt_vec;
        VectorEta_i_SimuOpt_T3 = d*(VectorK_i_SimuOpt_T3 + L -1) +1;
        DivisionSimuOpt_T3{1,index_r}(:,index1) = VectorK_i_SimuOpt_T3';


        %% To obtain the expectation of computational delay with/without cancellation

        SumDelayAPCC = 0;
        SumDelayAPCC_SimuOpt_T2 = 0;
        SumDelayAPCC_Cancellation = 0;
        SumDelayAPCC_Cancellation_SimuOpt_T3 = 0;
        for index3 = 1: RepeatTimes
            u1 = rand(r,N);
            VectorT = alpha0 / K - log(1-u1) / (mu0 * K);
            VectorTm = VectorT;

            DelayAPCCTemp = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i_MVD_T2); % to obtain the delay of the complete task
            SumDelayAPCC = SumDelayAPCC + DelayAPCCTemp;
            if sum(abs(VectorEta_i_SimuOpt_T2 - VectorEta_i_MVD_T2)) ~= 0 && r == r_BF
                DelayAPCCTemp2 = ComputeEndDelay_woCancel(VectorTm,r,VectorEta_i_SimuOpt_T2);
                SumDelayAPCC_SimuOpt_T2 = SumDelayAPCC_SimuOpt_T2 + DelayAPCCTemp2;
            end

            DelayAPCCTemp3 = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i_MVD_T3);
            SumDelayAPCC_Cancellation = SumDelayAPCC_Cancellation + DelayAPCCTemp3;
            if sum(abs(VectorEta_i_SimuOpt_T3 - VectorEta_i_MVD_T3)) ~= 0 && r == r_BF
                DelayAPCCTemp4 = ComputeEndDelay_wCancel(VectorTm,r,VectorEta_i_SimuOpt_T3);
                SumDelayAPCC_Cancellation_SimuOpt_T3 = SumDelayAPCC_Cancellation_SimuOpt_T3 + DelayAPCCTemp4;
            end
        end

        DelayAPCC = SumDelayAPCC / RepeatTimes;
        if sum(abs(VectorEta_i_SimuOpt_T2 - VectorEta_i_MVD_T2)) ~= 0 && r == r_BF
            DelayAPCC_SimuOpt_T2 = SumDelayAPCC_SimuOpt_T2 / RepeatTimes;
        else
            DelayAPCC_SimuOpt_T2 = DelayAPCC;
        end

        DelayAPCC_Cancellation = SumDelayAPCC_Cancellation / RepeatTimes;
        if sum(abs(VectorEta_i_SimuOpt_T3 - VectorEta_i_MVD_T3)) ~= 0 && r == r_BF
            DelayAPCC_SimuOpt_T3 = SumDelayAPCC_Cancellation_SimuOpt_T3 / RepeatTimes;
        else
            DelayAPCC_SimuOpt_T3 = DelayAPCC_Cancellation;
        end

        DelayAPCC_r_MVD_T2(index_r,index1) = DelayAPCC;
        DelayAPCC_r_SimuOpt_T2(index_r,index1) = DelayAPCC_SimuOpt_T2;
        DelayAPCC_r_MVD_T3(index_r,index1) = DelayAPCC_Cancellation;
        DelayAPCC_r_SimuOpt_T3(index_r,index1) = DelayAPCC_SimuOpt_T3;

    end % End of loop for KPrime

    Vec_r_DelayAPCCMinimum_MVD_T2(index_r) = min(DelayAPCC_r_MVD_T2(index_r,:));
    Vec_r_DelayAPCCMinimum_SimuOpt_T2(index_r) = min(DelayAPCC_r_SimuOpt_T2(index_r,:));
    Vec_r_DelayAPCCMinimum_MVD_T3(index_r) = min(DelayAPCC_r_MVD_T3(index_r,:));
    Vec_r_DelayAPCCMinimum_SimuOpt_T3(index_r) = min(DelayAPCC_r_SimuOpt_T3(index_r,:));

end% End of loop for r
disp(['Execution time of MVD_T2:',num2str(ExcuTimeMVD_T2)]);
disp(['Execution time of MVD_T3:',num2str(ExcuTimeMVD_T3)]);
disp(['Execution time of SimuOpt_T2:',num2str(ExcuTimeSimuOpt_T2)]);
disp(['Execution time of SimuOpt_T3:',num2str(ExcuTimeSimuOpt_T3)]);
[DelayAPCCMinimum_MVD_T2,temp2] = min(Vec_r_DelayAPCCMinimum_MVD_T2);
r_Optimal_MVD_T2 = Vector_r(temp2);
[DelayAPCCMinimum_MVD_T3,temp2] = min(Vec_r_DelayAPCCMinimum_MVD_T3);
r_Optimal_MVD_T3 = Vector_r(temp2);

%% Delay of Uncoded

% SumDelayUncoded = 0;
% for index3 = 1: RepeatTimes
%     u3 = rand(1,N);
%     VectorTPrime = alpha0 / N - log(1-u3) / (mu0 * N);
%     DelayUncodedTemp = max(VectorTPrime);
%     SumDelayUncoded = SumDelayUncoded + DelayUncodedTemp;
% end
% DelayUncoded = SumDelayUncoded / RepeatTimes *DelayUncoded;

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

end
DelayLCCMinimum = min(DelayLCC);

%% Delay of LCC-MMC

if L ==0
    for index_r = 1:length(Vector_r)
        r = Vector_r(index_r);
        for index1 = 1:length(VectorKPrime)

            KPrime = VectorKPrime(index1);
            K_LM = KPrime*r; % the number of division in LLC-MMC
            Eta_LM = d*(K_LM - 1) + 1; %L = 0
            SumDelayLCC_MMC = 0;
            for index3 = 1: RepeatTimes
                u3 = rand(r,N);
                VectorT_MMC = alpha0 / K_LM - log(1-u3) / (mu0 * K_LM);
                VectorT_MMC_sort = zeros(1,N*r);%把所有worker的延时便在同一行
                VectorT_MMC_sort(1:N) =  VectorT_MMC(1,:);
                % update the delay of each subtask for each worker
                for index4 = 2:r
                    VectorT_MMC(index4,:) = VectorT_MMC(index4,:) + VectorT_MMC(index4 - 1,:);
                    VectorT_MMC_sort(index4*N-N+1:index4*N) = VectorT_MMC(index4,:);
                end
                [SortedVectorT_MMC,~] = sort(VectorT_MMC_sort);
                DelayLCC_MMCTemp = SortedVectorT_MMC(Eta_LM);
                SumDelayLCC_MMC = SumDelayLCC_MMC + DelayLCC_MMCTemp;
            end
            DelayLCC_MMC(index_r,index1) = SumDelayLCC_MMC / RepeatTimes;

        end
        Vec_r_DelayLCC_MMC_Minimum(index_r) = min(DelayLCC_MMC(index_r,:));
    end
end

%% Simulation results

if Mark_simulation == 1

    Disp_r1 = r_BF;
    Disp_r2 = r_disp;
    DelayAPCC1 = DelayAPCC_r_MVD_T2(Vector_r == Disp_r1,:);
    DelayAPCC2 = DelayAPCC_r_MVD_T3(Vector_r == Disp_r1,:);

    if Mark_BF == 1
        DelayAPCC3 = DelayAPCC_r_SimuOpt_T2(Vector_r == r_BF,:);
        DelayAPCC4 = DelayAPCC_r_SimuOpt_T3(Vector_r == r_BF,:);
    end

    DelayAPCC5 = DelayAPCC_r_MVD_T2(Vector_r == Disp_r2,:);
    DelayAPCC6 = DelayAPCC_r_MVD_T3(Vector_r == Disp_r2,:);


    figure(1),hold on
    set(gcf,'Units','centimeter','Position',[5 5 20 10]);
    set(gca,'position',[0.08,0.14,0.90,0.84]);

    ytemp = 1000*DelayAPCC1(DelayAPCC1 ~= 0);
    xtemp = VectorKPrime(DelayAPCC1 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/o Cancel, r=' num2str(Disp_r1)])

    ytemp = 1000*DelayAPCC2(DelayAPCC2 ~= 0);
    xtemp = VectorKPrime(DelayAPCC2 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/ Cancel, r=' num2str(Disp_r1)])

    if Mark_BF == 1
        ytemp = 1000*DelayAPCC3(DelayAPCC3 ~= 0);
        xtemp = VectorKPrime(DelayAPCC3 ~= 0);
        plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['APCC w/o Cancel, BF, r=' num2str(r_BF)])

        ytemp = 1000*DelayAPCC4(DelayAPCC4 ~= 0);
        xtemp = VectorKPrime(DelayAPCC4 ~= 0);
        plot(xtemp,ytemp,'--','LineWidth',1.25,'DisplayName',['APCC w/ Cancel, BF, r=' num2str(r_BF)])
    end

    ytemp = 1000*DelayAPCC5(DelayAPCC5 ~= 0);
    xtemp = VectorKPrime(DelayAPCC5 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/o Cancel, r=' num2str(Disp_r2)])

    ytemp = 1000*DelayAPCC6(DelayAPCC6 ~= 0);
    xtemp = VectorKPrime(DelayAPCC6 ~= 0);
    plot(xtemp,ytemp,'LineWidth',1.25,'DisplayName',['APCC w/ Cancel, r=' num2str(Disp_r2)])

    % ytemp = DelayUncoded;
    % xtemp = VectorKPrime;
    % plot(xtemp,ytemp,'LineWidth',1.5,'DisplayName','Uncoded')

    ytemp = 1000*DelayLCC;
    xtemp = VectorKPrime;
    plot(xtemp,ytemp,'-','LineWidth',1.25,'DisplayName','LCC')

    if L==0
        ytemp = 1000*DelayLCC_MMC(Vector_r == Disp_r1,:);
        xtemp = VectorKPrime;
        plot(xtemp,ytemp,'-','LineWidth',1.25,'DisplayName',['LCC-MMC, r=' num2str(Disp_r1)])

        ytemp = 1000*DelayLCC_MMC(Vector_r == Disp_r2,:);
        xtemp = VectorKPrime;
        plot(xtemp,ytemp,'-','LineWidth',1.25,'DisplayName',['LCC-MMC, r=' num2str(Disp_r2)])
    end

    %     line([KPrimeMax+1 KPrimeMax+1],[0 max(DelayLCC)*1100],'linestyle','--','LineWidth', 1.25,'DisplayName','Maximum number of task segmentation, Lemma 1');

    %标签和坐标轴字体大小
    legend('FontSize',12,'Position',[0.2 0.7 0.3 0.3])
    set(gca, 'FontSize', 12);

    %     title('APCC w/o cancellation vs LCC')
    xlabel('K''=K/r, the number of task divisions')
    ylabel('Task completion delay (ms)')
    xlim([0 KPrimeMax+1])
    % ylim([0 800])

    %% plot the Histogram

    if Mark_Histogram == 1

        if L>0
            DelayList1 = NaN.*rand(2,7);
            %     DelayList1 = [rand(2,5);NaN.*ones(2,5)];

            %     LabelList = categorical({'N=100,L=10,d=2,r=4','N=100,L=10,d=2,r=16','N=200,L=20,d=4,r=4','N=200,L=20,d=4,r=16'});
            %     LabelList = reordercats(LabelList,{'N=100,L=10,d=2,r=4','N=100,L=10,d=2,r=16','N=200,L=20,d=4,r=4','N=200,L=20,d=4,r=16'});
            LabelList = categorical({'N=100,L=10,d=2','N=200,L=20,d=4'});
            LabelList = reordercats(LabelList,{'N=100,L=10,d=2','N=200,L=20,d=4'});

            if N==100 && L ==10 && d==2
                DelayList1(1,1) = 1000*DelayLCCMinimum;
                DelayList1(1,2) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_disp);
                DelayList1(1,3) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_disp);
                DelayList1(1,4) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_BF);
                DelayList1(1,5) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_BF);
                DelayList1(1,6) = 1000*Vec_r_DelayAPCCMinimum_SimuOpt_T2(Vector_r == r_BF);
                DelayList1(1,7) = 1000*Vec_r_DelayAPCCMinimum_SimuOpt_T3(Vector_r == r_BF);
            end

            %     if N==100 && L ==10 && d==2
            %         DelayList2(2,1) = DelayLCCMinimum;
            %         DelayList2(2,2) = Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == 16);
            %         DelayList2(2,3) = Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == 16);
            %     end

            if N==200 && L ==20 && d==4
                DelayList1(2,1) = 1000*DelayLCCMinimum;
                DelayList1(2,2) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_disp);
                DelayList1(2,3) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_disp);
                DelayList1(2,4) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_BF);
                DelayList1(2,5) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_BF);
                DelayList1(2,6) = 1000*Vec_r_DelayAPCCMinimum_SimuOpt_T2(Vector_r == r_BF);
                DelayList1(2,7) = 1000*Vec_r_DelayAPCCMinimum_SimuOpt_T3(Vector_r == r_BF);
            end

            %     if N==200 && L ==20 && d==4
            %         DelayList2(4,1) = DelayLCCMinimum;
            %         DelayList2(4,2) = Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == 16);
            %         DelayList2(4,3) = Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == 16);
            %     end

            figure(2),hold on
            set(gcf,'Units','centimeter','Position',[5 5 20 10]);
            set(gca,'position',[0.08,0.08,0.84,0.84]);
            b1 = bar(LabelList,DelayList1,0.7,'FaceColor','flat');

            for i = 1:size(DelayList1,2)
                xtips = b1(i).XEndPoints;
                ytips = b1(i).YEndPoints;
                labels = string(round(b1(i).YData,1));
                text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
            end
            c = hsv(21);
            for i = 1:size(DelayList1,2)
                set(b1(i),'FaceColor',c(2*i,:))
            end

            legend('LCC', ['APCC w/o Cancel, MVD, r=',num2str(r_disp)], ['APCC w/ Cancel, MVD, r=',num2str(r_disp)], ...,
                ['APCC w/o Cancel, MVD, r=',num2str(r_BF)], ['APCC w/ Cancel, MVD, r=',num2str(r_BF)],...,
                ['APCC w/o Cancel, Brute-Force, r=',num2str(r_BF)],['APCC w/ Cancel, Brute-Force, r=',num2str(r_BF)])
            
            %标签和坐标轴字体大小
            legend('FontSize',12,'Position',[0.2 0.7 0.3 0.3])
            set(gca, 'FontSize', 12);
            
            %     ylim([0 600])
            xlim([LabelList(1),LabelList(end)])
            ylabel('Minimum task completion delay (ms)')

        elseif L ==0

            DelayList1 = NaN.*rand(2,7);

            LabelList = categorical({'N=50,L=0,d=2','N=100,L=0,d=4'});
            LabelList = reordercats(LabelList,{'N=50,L=0,d=2','N=100,L=0,d=4'});

            if N==100 && L ==0 && d==2
                DelayList1(1,1) = 1000*DelayLCCMinimum;
                DelayList1(1,2) = 1000*Vec_r_DelayLCC_MMC_Minimum(Vector_r==r_BF);
                DelayList1(1,3) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_BF);
                DelayList1(1,4) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_BF);
                DelayList1(1,5) = 1000*Vec_r_DelayLCC_MMC_Minimum(Vector_r==r_disp);
                DelayList1(1,6) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_disp);
                DelayList1(1,7) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_disp);
            end

            if N==200 && L ==0 && d==4
                load Delay_dat.mat % N==100 && L ==0 && d==4时读取之前存的数据
                DelayList1(2,1) = 1000*DelayLCCMinimum;
                DelayList1(2,2) = 1000*Vec_r_DelayLCC_MMC_Minimum(Vector_r==r_BF);
                DelayList1(2,3) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_BF);
                DelayList1(2,4) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_BF);
                DelayList1(2,5) = 1000*Vec_r_DelayLCC_MMC_Minimum(Vector_r==r_disp);
                DelayList1(2,6) = 1000*Vec_r_DelayAPCCMinimum_MVD_T2(Vector_r == r_disp);
                DelayList1(2,7) = 1000*Vec_r_DelayAPCCMinimum_MVD_T3(Vector_r == r_disp);
            end
            save Delay_dat.mat DelayList1

            figure(2),hold on
            grid on
            set(gcf,'Units','centimeter','Position',[5 5 20 10]);
            set(gca,'position',[0.08,0.08,0.91,0.90]);
            b1 = bar(LabelList,DelayList1,0.7,'FaceColor','flat');

            for i = 1:size(DelayList1,2)
                xtips = b1(i).XEndPoints;
                ytips = b1(i).YEndPoints;
                labels = string(round(b1(i).YData,1));
                text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom')
            end
            c = hsv(2*size(DelayList1,2)+1);
            for i = 1:size(DelayList1,2)
                set(b1(i),'FaceColor',c(2*i-1,:))
            end

            legend('LCC', ['LCC-MMC, r=' num2str(r_BF)],['APCC w/o Cancel, MVD, r=',num2str(r_BF)], ['APCC w/ Cancel, MVD, r=',num2str(r_BF)], ...,
                ['LCC-MMC, r=' num2str(r_disp)], ['APCC w/o Cancel, MVD, r=',num2str(r_disp)], ['APCC w/ Cancel, MVD, r=',num2str(r_disp)])

            %标签和坐标轴字体大小
            legend('FontSize',12,'Position',[0.2 0.7 0.3 0.3])
            set(gca, 'FontSize', 12);

            ylim([0 500])
            xlim([LabelList(1),LabelList(end)])
            ylabel('Minimum task completion delay (ms)')
        end
    end

    if Mark_rL == 1
        figure(3),hold on
        set(gcf,'Units','centimeter','Position',[5 5 20 10]);
        set(gca,'position',[0.08,0.08,0.84,0.84]);

        plot(Vector_r,1000.*Vec_r_DelayAPCCMinimum_MVD_T2,'-s','LineWidth',1.25,'DisplayName',['N=',num2str(N),', L=',num2str(L),', d=',num2str(d), ', w/o Cancel'])
        plot(Vector_r,1000.*Vec_r_DelayAPCCMinimum_MVD_T3,'-s','LineWidth',1.25,'DisplayName',['N=',num2str(N),', L=',num2str(L),', d=',num2str(d), ', w/ Cancel'])

        %     title('APCC w/o cancellation vs LCC')
        xlabel('r, the hierarchical division number of sets')
        ylabel('Minimum Average Delay (ms)')
        %     xlim([0 KPrimeMax+1])
        %     ylim([0 800])
    end

end

toc(timerval)