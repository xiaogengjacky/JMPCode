clear all
close all

%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%

%% Basic Setup

%addpath /home/jy489/Matlab;

D_Constant;

Dpara = 10;  %0 - nonparametric probability weighting; 
            %1 - parametric probability weighting a+b*p; 
            %10 - non parametric with heterogeniety in loss aversion
            %11 - parametric a+b*p with heterogeniety in loss aversion
DMug = 0;  %0 - Money: one dimension; 1 - Mug: tow dimensions
Fform = 0;  %0 - 1 dimension; 1 - 2 dimensions
ErrDist = 0;  %0 - additive errors;  1 - multiplicative errors;
chain_length = 1e5;

load WTAWTP.mat;

if DMug == 0
    data = Money;
elseif DMug == 1
    data = Mug;
end

PT = [1 1 1 1 1 1 1 1 1];  %Specify which probabilities to include
TT = [1 1 1 1];  %Specify which treatments to include

if Dpara == 0 || Dpara == 10
    NV = 12 + floor(Dpara/10); % number of parameters
elseif Dpara == 1 || Dpara == 11
    NV = 5 + floor(Dpara/10); % number of parameters
end

%data = ManipulateObv(data,PT,TT);
N = size(data,1);
flag = [Dpara Fform ErrDist];
%% intial guess
rng(12345);

coeffs_mcmc    =   zeros(NV,chain_length);
LikMag         =   zeros(1,chain_length);
update         =   zeros(1,chain_length);


if Dpara == 0
%     coeffs_mcmc(:,1) = DMug*[1.5; 0; 7.2; 0.09; 0.16; 0.22; 0.39; 0.52; 0.61; 0.70; 0.73; 0.75]...
%         +(1-DMug)*[2.5; 0; 7.8; 0.05; 0.14; 0.19; 0.31; 0.47; 0.58; 0.71; 0.74; 0.76];
    coeffs_mcmc(:,1) = DMug*[1; 0; 10; 0; 0.005; 0.025; 0.125; 0.25; 0.375; 0.475; 0.495; 0.5]...
         +(1-DMug)*[1; 0; 5; 0; 0.01; 0.05; 0.25; 0.5; 0.75; 0.95; 0.99; 1];

elseif Dpara == 1
    coeffs_mcmc(:,1) = DMug*[1.5; 0; 7.3; 0.17; 0.58]...
        +(1-DMug)*[2.5; 0; 7.9; 0.12; 0.63];
elseif Dpara == 10
%     coeffs_mcmc(:,1) = DMug*[1.5; 0; 7.2; 0.09; 0.16; 0.22; 0.39; 0.52; 0.61; 0.70; 0.73; 0.75; 0]...
%         +(1-DMug)*[2.5; 0; 7.8; 0.05; 0.14; 0.19; 0.31; 0.47; 0.58; 0.71; 0.74; 0.76; 0];
    coeffs_mcmc(:,1) = DMug*[1; 0; 10; 0; 0.005; 0.025; 0.125; 0.25; 0.375; 0.475; 0.495; 0.5; 0]...
         +(1-DMug)*[1; 0; 5; 0; 0.01; 0.05; 0.25; 0.5; 0.75; 0.95; 0.99; 1; 0];
elseif Dpara == 11
    coeffs_mcmc(:,1) = DMug*[1.5; 0; 7.3; 0.17; 0.58; 0]...
        +(1-DMug)*[2.5; 0; 7.; 0.12; 0.63; 0];
end

LikMag(1,1) = 0;

if floor(Dpara/10)==1
    hermite_rule ( 20, 0, 0.5, 'hermite_for_mcmc' );
    ws = dlmread('hermite_for_mcmc_w.txt');
    ns = dlmread('hermite_for_mcmc_x.txt');
    LL_incumbent = LLike_Integrate(coeffs_mcmc(:,1), data, ws, ns, flag, PT);
elseif floor(Dpara/10)==0
    LL_incumbent = LLike(coeffs_mcmc(:,1), data, flag, PT);
end

if Fform == 0
    step = [0.1; 0.1; 0.1];
     step = [step; 0.01*ones(9,1)];
%Constraint on probability weighting
%    step = [step; 0*ones(9,1)];
elseif Fform == 1
    step = [0.1; 0.1; 0.2];
    step = [step; 0.005*ones(9,1)];
end

if Dpara >=10 
    step = [step; 0.1];
end

%% Update
tic
for ii=2:chain_length
 
    if (ii/10000-round(ii/10000))==0
        disp(ii)
        toc
    end
    
    %k = mod(ii,NV)+1;
    %k = randi(NV,1);
    k = randi(2,[NV, 1]) - ones(NV, 1);
    %B(k) = step(k);
    B = k.*step;
    
    coeff_incumbent = coeffs_mcmc(:,ii-1);
    coeff_candidate = coeffs_mcmc(:,ii-1)+B.*randn(NV,1);
    
    if update(ii-1)==1
        if floor(Dpara/10)==1
            LL_incumbent  = LLike_Integrate(coeff_incumbent, data, ws, ns, flag, PT);
        else            
            LL_incumbent  = LLike(coeff_incumbent, data, flag, PT);
        end
    end
    
    if floor(Dpara/10)==1
        LL_candidate  = LLike_Integrate(coeff_candidate, data, ws, ns, flag, PT);
    else
        LL_candidate  = LLike(coeff_candidate, data, flag, PT);
    end
        
    criterion = exp(-LL_candidate+LL_incumbent);

    if criterion>rand(1)
        coeffs_mcmc(:,ii)= coeff_candidate;
        LikMag(:,ii) = LL_candidate;
        update(ii) = 1;
    else      
        coeffs_mcmc(:,ii)= coeff_incumbent;
        LikMag(:,ii) = LL_incumbent;
        update(ii) = 0;
    end
    
end

coeffs_mcmc(2,:) = exp(coeffs_mcmc(2,:)); %convert log(sigma) to sigma

if floor(Dpara/10)==1
    coeffs_mcmc(size(coeffs_mcmc,1),:) = exp(coeffs_mcmc(size(coeffs_mcmc,1),:));
end

save ('mcmc.mat');
%% Diagnostics

if Fform == 0
    figure;

    subplot 221
    plot(1:chain_length,coeffs_mcmc(1,:),'r','linewidth',2)
    axis([0 chain_length 0 5]);
    title('Loss Aversion Parameter (\beta)', 'fontsize', 12)

    subplot 222
    plot(1:chain_length,coeffs_mcmc(2,:),'r','linewidth',2)
    axis([0 chain_length 0 5]);
    title('Std. Dev of Error Term in Econometric Model (\sigma)', 'fontsize', 12)

    subplot 223
    plot(1:chain_length,coeffs_mcmc(3,:),'r','linewidth',2)
    axis([0 chain_length 0 12]);
    title('Mean Value (X)', 'fontsize', 12)

%     subplot 224
%     plot(1:chain_length,coeffs_mcmc(13,:),'r','linewidth',2)
%     axis([0 chain_length 0 12]);
%     title('Standard Deviation of Value (X)', 'fontsize', 12)

    figure;
    subplot 331
    plot(1:chain_length,coeffs_mcmc(4,:),'r','linewidth',2)
    axis([0 chain_length 0 1]);
    title('Decision Weight \pi(0)')

    subplot 332
    plot(1:chain_length,coeffs_mcmc(5,:),'r','linewidth',2)
    axis([0 chain_length 0 1]);
    title('Decision Weight \pi(0.01)')

    if Dpara == 0 || Dpara == 10
        subplot 333
        plot(1:chain_length,coeffs_mcmc(6,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.05)')

        subplot 334
        plot(1:chain_length,coeffs_mcmc(7,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.25)')

        subplot 335
        plot(1:chain_length,coeffs_mcmc(8,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.50)')

        subplot 336
        plot(1:chain_length,coeffs_mcmc(9,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.75)')

        subplot 337
        plot(1:chain_length,coeffs_mcmc(10,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.95)')

        subplot 338
        plot(1:chain_length,coeffs_mcmc(11,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(0.99)')

        subplot 339
        plot(1:chain_length,coeffs_mcmc(12,:),'r','linewidth',2)
        axis([0 chain_length 0 1]);
        title('Decision Weight \pi(1)')

    end
    
elseif Fform ==1 

    figure;

    subplot 221
    plot(1:chain_length,coeffs_mcmc(1,:),'b','linewidth',2)
    axis([0 chain_length 0 5]);
    title('Loss Aversion Parameter (\beta)', 'fontsize', 14)

    subplot 222
    plot(1:chain_length,coeffs_mcmc(2,:),'b','linewidth',2)
    axis([0 chain_length 0 5]);
    title('Std. Dev of Error Term in Econometric Model (\sigma)', 'fontsize', 14)

    subplot 223
    plot(1:chain_length,coeffs_mcmc(3,:),'b','linewidth',2)
    axis([0 chain_length 0 30]);
    title('m/a', 'fontsize', 12)

    subplot 224
    plot(1:chain_length,coeffs_mcmc(13,:),'b','linewidth',2)
    axis([0 chain_length 0 12]);
    title('Standard Deviation of m/a', 'fontsize', 14)

    figure;
    subplot 331
    plot(1:chain_length,coeffs_mcmc(4,:),'b','linewidth',2)
    axis([0 chain_length 0 10]);
    title('a\pi(0)', 'fontsize', 14)

    subplot 332
    plot(1:chain_length,coeffs_mcmc(5,:),'b','linewidth',2)
    axis([0 chain_length 0 10]);
    title('a\pi(0.01)', 'fontsize', 14)

    if Dpara == 0 || Dpara == 10
        subplot 333
        plot(1:chain_length,coeffs_mcmc(6,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.05)', 'fontsize', 14)

        subplot 334
        plot(1:chain_length,coeffs_mcmc(7,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.25)', 'fontsize', 14)

        subplot 335
        plot(1:chain_length,coeffs_mcmc(8,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.50)', 'fontsize', 14)

        subplot 336
        plot(1:chain_length,coeffs_mcmc(9,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.75)', 'fontsize', 14)

        subplot 337
        plot(1:chain_length,coeffs_mcmc(10,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.95)', 'fontsize', 14)

        subplot 338
        plot(1:chain_length,coeffs_mcmc(11,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(0.99)', 'fontsize', 14)

        subplot 339
        plot(1:chain_length,coeffs_mcmc(12,:),'b','linewidth',2)
        axis([0 chain_length 0 1]);
        title('a\pi(1)', 'fontsize', 14)

    end
    
end


figure;

plot(1:chain_length,cumsum(update(1,:))./(1:chain_length),'r','linewidth',2)
axis([0 chain_length 0 1])
title('Average Occurence of Updating', 'fontsize', 14)
xlabel('iteration')

%% Results Processing

mcmc = coeffs_mcmc(:,chain_length/4:end);
mcmc = mcmc(:,1:200:end);

Rpoint = 1; %Reference point - between 0 and 1
if Fform == 1
    mcmc(3,:) = mcmc(3,:).*Rpoint;
    for jj = 4:12
        mcmc(jj,:) = mcmc(jj,:)./Rpoint;
    end
    
    if Dpara == 10
        mcmc(13,:) = mcmc(13,:).*Rpoint;
    end
end

if Dpara == 0
    Coe_els = mean(mcmc(1:12,:),2)
    Coe_els_std = std(mcmc(1:12,:),0,2)
elseif Dpara == 10
    Coe_els = mean(mcmc(1:13,:),2)
    Coe_els_std = std(mcmc(1:13,:),0,2)
elseif Dpara == 1
    Coe_els = mean(mcmc(1:5,:),2)
    Coe_els_std = std(mcmc(1:5,:),0,2)
elseif Dpara == 11
    Coe_els = mean(mcmc(1:6,:),2)
    Coe_els_std = std(mcmc(1:6,:),0,2)
end

[wtpg, wtag, wtpl, wtal, pi_p, MSE, RSQ] = RDisp(mcmc,data,flag,Rpoint)

p = [0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];

figure;
hold off;
plot(p,wtpg,'bs-', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'w');
axis([0 1 0 8]);
title('Model Prediction of Money($5) Experiment', 'fontsize', 14);
xlabel('Objective Probability', 'fontsize', 12);
ylabel('Valuation', 'fontsize', 12);
hold on;
plot(p,wtag,'b^:', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'w');
plot(p,wtpl,'ks-.', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','g', 'MarkerEdgeColor', 'w');
plot(p,wtal,'k^--', 'LineWidth',1.5,'MarkerSize',8,  'MarkerFaceColor','g', 'MarkerEdgeColor', 'w');
legend('WTP-Gain','WTA-Gain','WTP-Loss','WTA-Loss','Location','SouthEast');
hold off;

figure;
plot(p, pi_p, 'bo-', 'LineWidth', 1.5, 'MarkerSize',7);
axis([0 1 0 1]);
title('Nonparametric Probability Weighting Function', 'fontsize', 14);
xlabel('Objective Probability', 'fontsize', 12);
ylabel('Decision Weight','fontsize',12);
hold on;
plot(p, p, 'k:');
legend('Weighted Probability','45 Degree Line','Location','SouthEast');

figure;
plot(p, pi_p./0.4, 'ko-.', 'LineWidth', 1.5, 'MarkerSize',4);
axis([0 1 0 1]);
title('Nonparametric Probability Weighting Function', 'fontsize', 14);
xlabel('Objective Probability', 'fontsize', 12);
ylabel('Decision Weight','fontsize',12);
hold on
%plot(p, pi_p./0.4, 'ro-', 'LineWidth', 1.5, 'MarkerSize',4);
plot(p, pi_p./0.5, 'ro:', 'LineWidth', 1.5, 'MarkerSize',4);
plot(p, pi_p./0.6, 'bo-', 'LineWidth', 1.5, 'MarkerSize',4);
plot(p, pi_p./0.7, 'mo--', 'LineWidth', 1.5, 'MarkerSize',4);
plot(p, p, 'k:');
legend('Weighted Probability (a=0.4)','Weighted Probability (a=0.5)',...
    'Weighted Probability (a=0.6)','Weighted Probability (a=0.7)','45 Degree Line','Location','SouthEast');