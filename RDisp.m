function [x1, x2, x3, x4, pi_p, MSE, RS, LL]=RDisp(para,data,flag,PT)
D_Constant;
base = 4;
Sig = para(2);

L = mean(para(3,:),2);
L1 = L;
L2 = L;
L3 = L;
L4 = L;

Dpara = flag(1);
Dform = flag(2);
a = PT;

gam = mean(para(base,:),2);
if Dpara<=1 || Dpara == 11
    del = mean(para(base+1,:),2);
end

onesp = ones(1,9);
onespp = ones(size(data,1),9);

if Dpara == 0 || Dpara == 10
    pi_p = mean(para(base:base+8,:),2);
    %pi_p(PT'==0) = [];
elseif Dpara == 1 || Dpara == 11
    pi_p = gam*onesp + del*p;
end
pi_p = pi_p';
h = 1./pi_p-onesp;
g = mean(para(1,:),2);

if Dform == 0
    x1  =  (L1./(onesp+h.*g))';
    x2  =  (pi_p*L2)';
    x3  =  (pi_p*L3)';
    x4  =  (L4./(onesp+h./g))';  
elseif Dform == 1
    x1  =  (pi_p*L1./g)';
    x2  =  (pi_p*L2+pi_p.^2*(g-1)*L2*a)';
    x3  =  (pi_p*L3+pi_p.^2*(1/g-1)*L3*a)';
    x4  =  (pi_p*L4.*g)';
end

pi_p2 = repmat(pi_p, size(data,1),1);

if Dform == 0
    h = 1./pi_p2-onespp;
    Prediction = (L1./(onespp+h*g)).*repmat(data(:,T1),1,9) + (pi_p2*L2).*repmat(data(:,T2),1,9) +...
        (pi_p2*L3).*repmat(data(:,T3),1,9) + (L4./(onespp+h/g)).*repmat(data(:,T4),1,9);
elseif Dform == 1
    Prediction = (pi_p2*L1./g).*repmat(data(:,T1),1,9) + (pi_p2*L2+pi_p2.^2*(g-1)*L2*a).*repmat(data(:,T2),1,9) +...
        (pi_p2*L3+pi_p2.^2*(1/g-1)*L3*a).*repmat(data(:,T3),1,9) + (pi_p2*L4.*g).*repmat(data(:,T4),1,9);
end

RSS = (data(:,BID1:BID9) - Prediction);
LL = normpdf(RSS, 0, Sig);
LL = log(LL);
LL = sum(sum(LL,2),1);
RSS = (data(:,BID1:BID9) - Prediction).^2;
RSS = sum(RSS(:));
TSS = data(:,BID1:BID9);
TSS = TSS(:);
TSS = var(TSS)*size(data(:,BID1:BID9),1)*size(data(:, BID1:BID9),2);
%TSS = sum(TSS(:));
RS = 1-RSS./TSS;
MSE = RSS/(size(data(:,BID1:BID9),1)*size(data(:, BID1:BID9),2));
