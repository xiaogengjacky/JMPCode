clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program jointly estimate all decision weights and variance of bids
%% Basic Setup
hermite_rule ( 20, 0, 0.5, 'Hermite_for_7120' );
ws =   dlmread('hermite_for_7120_w.txt');
ns   =   dlmread('hermite_for_7120_x.txt');

D_Constant;

% parameters setup
Model_s = 0; %0 - EV; 1 - EU; 2 - RDU; 3 - RDU with PW
Dpara = 0;  %0 - nonparametric probability weighting; 
            %1 - parametric probability weighting a+b*p; 
            %10 - non parametric with heterogeniety in loss aversion
            %11 - parametric a+b*p with heterogeniety in loss aversion
DMug = 1;  %0 - Money: one dimension; 1 - Mug: tow dimensions
Fform = 1;  %0 - 1 dimension; 1 - 2 dimensions
ErrDist = 0;  %0 - additive errors;  1 - multiplicative errors;
Drestrict = 0;  %0 - no restriction;  1 - imposing the restriction that pi(0)=0 and pi(1)=1
Nboot = 1;  % number of bootstrap iterations; 1 - no bootstrap 

Rpoint = 0.5;

options = optimset('Algorithm', 'interior-point', 'Display','notify',... 
            'TolX',1e-6, 'TolFun', 1e-6,... 
            'MaxFunEvals',1e5,'MaxIter',2e3);
load WTAWTP.mat;

if DMug == 0
    data = Money_C;
    data_v = Money_V;
elseif DMug == 1
    data = Mug_C;
    data_v = Mug_V;
end

PT = [1 1 1 1 1 1 1 1 1];  %Specify which probabilities to include
TT = [1 1 1 1];  %Specify which treatments to include

%data = ManipulateObv(data,PT,TT);
N = size(data,1);
%Het_Sig;
    
%%  Estimation - MLE
rng(12345);
flag = [Dpara Fform ErrDist Rpoint];
%Prameter Vector - beta, sigma, ubar, probability weighting functions
% for hetergeniety specification, there is an extra parameter for hetero

if Dpara == 0 || Dpara == 10
    if Model_s == 3
        para0  =  [10*rand(1) 3 10*rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1) rand(1)];
        Lower  =  [0 -Inf 0 0 0 0 0 0 0 0 0 0];
        Upper  =  [50 Inf 12+12*Fform];
        Upper  =  [Upper PT];
    elseif Model_s == 2
%   Estimate restricted reference dependent utility model    
        para0  =  [10*rand(1) 3 10*rand(1) 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Lower  =  [0 -Inf 0 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Upper  =  [50 Inf 12+12*Fform];
        Upper  =  [Upper 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
    elseif Model_s == 1
%   Estimate expected utility model
        para0  =  [1 3 10*rand(1) 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Lower  =  [1 -Inf 0 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Upper  =  [1 Inf 12+12*Fform];
        Upper  =  [Upper 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
    elseif Model_s == 0
        para0  =  [1 3 5 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Lower  =  [1 -Inf 7.99 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
        Upper  =  [1 Inf 7.99];
        Upper  =  [Upper 0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];    
    end
        
    if Drestrict==1
        Upper(4) = 0;
        Lower(12) = 1;
    end
    
elseif Dpara == 1 || Dpara == 11
    para0  =  [10*rand(1) 10*rand(1) rand(1) rand(1) rand(1)];
    Lower  =  [0 -Inf 0 0 0];
    Upper  =  [20 Inf 12 2 2];
end

if Dpara >= 10
    para0  =  [para0 10*rand(1)];
    Lower  =  [Lower  0];
    Upper  =  [Upper  5];
end

parfor jj = 1:Nboot
    if(jj/10-round(jj/10))==0

    end
        
    weights = ws;
    nodes = ns;
    
    data_boot = data;
    if Nboot>1
        data_boot = data(randi(N,N,1),:);
    end
        
    if Dpara >= 10
        problem = createOptimProblem('fmincon', 'objective', @(para) LLike_Integrate(para,data_boot,weights,nodes,flag,PT), 'x0', para0, 'lb', Lower, 'ub', Upper, 'options', options);
    else
        problem = createOptimProblem('fmincon', 'objective', @(para) LLike(para,data_boot,flag,PT), 'x0', para0, 'lb', Lower, 'ub', Upper, 'options', options);
    end
    
    %MS = MultiStart('UseParallel', 'always', 'Display', 'final');
    %[x, f1, f2]=run(MS, problem,20);
    GS = GlobalSearch;
    [x, f1, f2] = run(GS, problem);
    if f2 >0 
        para_final(:,jj) = x;
    else
        para_final(:,jj) = NaN; %para_final(:,jj-1);
    end
    success(:,jj) = [f1; f2];
    
end

para_final(2,:) = exp(para_final(2,:)); %convert log(sigma) to sigma
save ('mle.mat');
%% Results Processing

coeff = para_final;

Rpoint = 0.5;

if Fform == 1
    coeff(3,:) = coeff(3,:).*Rpoint;
    for jj = 4:12
        coeff(jj,:) = coeff(jj,:)./Rpoint;
    end
    
    if Dpara == 10
        coeff(13,:) = coeff(13,:).*Rpoint;
    end
end

if Dpara == 0
    Coe_els = mean(coeff(1:12,:),2)
    Coe_els_std = std(coeff(1:12,:),0,2)
elseif Dpara == 10
    Coe_els = mean(coeff(1:13,:),2)
    Coe_els_std = std(coeff(1:13,:),0,2)
elseif Dpara == 1
    Coe_els = mean(coeff(1:5,:),2)
    Coe_els_std = std(coeff(1:5,:),0,2)
elseif Dpara == 11
    Coe_els = mean(coeff(1:6,:),2)
    Coe_els_std = std(coeff(1:6,:),0,2)
end

[wtpg, wtag, wtpl, wtal, pi_p, MSE, RSQ, LL] = RDisp(coeff,data_v,flag,Rpoint)

p = [0 0.01 0.05 0.25 0.5 0.75 0.95 0.99 1];
%pi_p = pi_p*(1/pi_p(9));

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

%% emial notification
% myaddress = 'jbyanjacky@gmail.com';
% mypassword = 'xiaogeng2363371';
% setpref('Internet','E_mail',myaddress);
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username',myaddress);
% setpref('Internet','SMTP_Password',mypassword);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', ...
%                   'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% sendmail(myaddress, 'Task update', 'Matlab task completed ', 'mle.mat');