function data = ManObv(fulldata, probability, treatment)
Def_Constant;
data = fulldata;
data((data(:,P00)==probability(1)+1),:) = [];
data((data(:,P01)==probability(2)+1),:) = [];
data((data(:,P05)==probability(3)+1),:) = [];
data((data(:,P25)==probability(4)+1),:) = [];
data((data(:,P50)==probability(5)+1),:) = [];
data((data(:,P75)==probability(6)+1),:) = [];
data((data(:,P95)==probability(7)+1),:) = [];
data((data(:,P99)==probability(8)+1),:) = [];
data((data(:,P100)==probability(9)+1),:) = [];

data((data(:,T1)==treatment(1)+1),:) = [];
data((data(:,T2)==treatment(2)+1),:) = [];
data((data(:,T3)==treatment(3)+1),:) = [];
data((data(:,T4)==treatment(4)+1),:) = [];

%data((data(:,BID)>10),:) = [];





