clear all
close all

%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%

% this code illustrates how mcmc - metropolis-hastings algorithm works for
% the simple case of estimating the Poisson claim rate and the mixture
% parameter.

%% Basic Setup

%addpath /home/jy489/Matlab;

Def_Constant;

Dpara = 10;  %0 - nonparametric probability weighting; 
            %1 - parametric probability weighting a+b*p; 
            %10 - non parametric with heterogeniety in loss aversion
            %11 - parametric a+b*p with heterogeniety in loss aversion
Mug = 0;  %0 - Money: one dimension; 1 - Mug: tow dimensions
Fform = 0;  %0 - 1 dimension; 1 - 2 dimensions
ErrDist = 0;  %0 - additive errors;  1 - multiplicative errors;
Rpoint = 0.5; %Reference point - between 0 and 1
chain_length = 2e4;

load MatlabData;

if Mug == 0
    data = MatlabMoney;
elseif Mug == 1
    data = MatlabMug;
end

PT = [1 1 1 1 1 1 1 1 1];  %Specify which probabilities to include
TT = [1 1 1 1];  %Specify which treatments to include

if Dpara == 0 || Dpara == 10
    NV = 12 + floor(Dpara/10); % number of parameters
elseif Dpara == 1 || Dpara == 11
    NV = 5 + floor(Dpara/10); % number of parameters
end

data = ManipulateObv(data,PT,TT);
N = size(data,1);
flag = [Dpara Fform ErrDist Rpoint];
%% intial guess
rng(12345);

coeffs_mcmc    =   zeros(NV,chain_length);
update         =   zeros(1,chain_length);


if Dpara == 0
    coeffs_mcmc(:,1) = Mug*[1.5; 0; 7.2; 0.09; 0.16; 0.22; 0.39; 0.52; 0.61; 0.70; 0.73; 0.75]...
        +(1-Mug)*[2.5; 0; 7.8; 0.05; 0.14; 0.19; 0.31; 0.47; 0.58; 0.71; 0.74; 0.76];
elseif Dpara == 1
    coeffs_mcmc(:,1) = Mug*[1.5; 0; 7.3; 0.17; 0.58]...
        +(1-Mug)*[2.5; 0; 7.9; 0.12; 0.63];
elseif Dpara == 10
    coeffs_mcmc(:,1) = Mug*[1.5; 0; 7.2; 0.09; 0.16; 0.22; 0.39; 0.52; 0.61; 0.70; 0.73; 0.75; 0]...
        +(1-Mug)*[2.5; 0; 7.8; 0.05; 0.14; 0.19; 0.31; 0.47; 0.58; 0.71; 0.74; 0.76; 0];
elseif Dpara == 11
    coeffs_mcmc(:,1) = Mug*[1.5; 0; 7.3; 0.17; 0.58; 0]...
        +(1-Mug)*[2.5; 0; 7.; 0.12; 0.63; 0];
end


if floor(Dpara/10)==1
    hermite_rule ( 20, 0, 0.5, 'hermite_for_mcmc' );
    ws = dlmread('hermite_for_mcmc_w.txt');
    ns = dlmread('hermite_for_mcmc_x.txt');
    LL_incumbent = LogLike_Integrate(coeffs_mcmc(:,1), data, ws, ns, flag);
elseif floor(Dpara/10)==0
    LL_incumbent = LogLike(coeffs_mcmc(:,1), data, flag);
end

step = [0.05; 0.1; 0.1];
step = [step; 0.01*ones(NV-3,1)];
%% Update
tic
for ii=2:chain_length
 
    if (ii/1000-round(ii/1000))==0
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
            LL_incumbent  = LogLike_Integrate(coeff_incumbent, data, ws, ns, flag);
        else            
            LL_incumbent  = LogLike(coeff_incumbent, data, flag);
        end
    end
    
    if floor(Dpara/10)==1
        LL_candidate  = LogLike_Integrate(coeff_candidate, data, ws, ns, flag);
    else
        LL_candidate  = LogLike(coeff_candidate, data, flag);
    end
        
    criterion = exp(-LL_candidate+LL_incumbent);

    if criterion>rand(1)
        coeffs_mcmc(:,ii)= coeff_candidate;
        update(ii) = 1;
    else      
        coeffs_mcmc(:,ii)= coeff_incumbent;
        update(ii) = 0;
    end
    
end

coeffs_mcmc(2,:) = exp(coeffs_mcmc(2,:)); %convert log(sigma) to sigma

if floor(Dpara/10)==1
    coeffs_mcmc(size(coeffs_mcmc,1),:) = exp(coeffs_mcmc(size(coeffs_mcmc,1),:));
end

save ('mcmc.mat');
%% Diagnostics

figure;

subplot 221
plot(1:chain_length,coeffs_mcmc(1,:),'b','linewidth',2)
axis([0 chain_length 0 5]);
title('\beta Markov chain')
%legend('MCMC','MLE','true')

subplot 222
plot(1:chain_length,coeffs_mcmc(2,:),'b','linewidth',2)
axis([0 chain_length 0 5]);
title('\sigma Markov chain')

subplot 223
plot(1:chain_length,coeffs_mcmc(3,:),'b','linewidth',2)
axis([0 chain_length 0 12]);
title('x Markov chain')

% subplot 224
% plot(1:chain_length,coeffs_mcmc(13,:),'b','linewidth',2)
% title('\sigma v Markov chain')

figure;
subplot 331
plot(1:chain_length,coeffs_mcmc(4,:),'b','linewidth',2)
axis([0 chain_length 0 1]);
title('\pi(0) Markov chain')

subplot 332
plot(1:chain_length,coeffs_mcmc(5,:),'b','linewidth',2)
axis([0 chain_length 0 1]);
title('\pi(0.01) Markov chain')

if Dpara == 0 || Dpara == 10
    subplot 333
    plot(1:chain_length,coeffs_mcmc(6,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.05) Markov chain')

    subplot 334
    plot(1:chain_length,coeffs_mcmc(7,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.25) Markov chain')

    subplot 335
    plot(1:chain_length,coeffs_mcmc(8,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.50) Markov chain')

    subplot 336
    plot(1:chain_length,coeffs_mcmc(9,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.75) Markov chain')

    subplot 337
    plot(1:chain_length,coeffs_mcmc(10,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.95) Markov chain')

    subplot 338
    plot(1:chain_length,coeffs_mcmc(11,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(0.99) Markov chain')

    subplot 339
    plot(1:chain_length,coeffs_mcmc(12,:),'b','linewidth',2)
    axis([0 chain_length 0 1]);
    title('\pi(1) Markov chain')

end

figure;

plot(1:chain_length,cumsum(update(1,:))./(1:chain_length),'linewidth',2)
axis([0 chain_length 0 1])
title('Average occurence of updating of \beta')
xlabel('iteration')

%% Results Processing
mcmc = coeffs_mcmc(:,chain_length/4:chain_length);
mcmc = mcmc(:,1:200:end);

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

[wtpg, wtag, wtpl, wtal, p, pi_p, MSE, RSQ] = ResDisp(mcmc,data,flag,PT)

figure;
hold off;
plot(p,wtpg,'bs-', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'w');
axis([0 1 0 8]);
title('Model Prediction of Money($5) Experiment', 'fontsize', 12);
xlabel('Probability');
ylabel('Valuation');
hold on;
plot(p,wtag,'b^:', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'w');
plot(p,wtpl,'ks-.', 'LineWidth',1.5,'MarkerSize',8, 'MarkerFaceColor','g', 'MarkerEdgeColor', 'w');
plot(p,wtal,'k^--', 'LineWidth',1.5,'MarkerSize',8,  'MarkerFaceColor','g', 'MarkerEdgeColor', 'w');
legend('WTP-Gain','WTA-Gain','WTP-Loss','WTA-Loss','Location','SouthEast');
hold off;

figure;
plot(p, pi_p, 'bo-', 'LineWidth', 1.5, 'MarkerSize',7);
title('Nonparametric Probability Weighting Function', 'fontsize', 12);
hold on;
plot(p, p, 'k:');
legend('Weighted Probability','45 Degree Line','Location','SouthEast');

