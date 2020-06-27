function [U_new,center,Kdist] = mystepfcm(data, cluster_n, expo,center)

epsilon=0.00000001;

Kdist=myDistance(center,cluster_n,data);
dist=1-Kdist;
dist=dist+epsilon;

tmp = dist.^(-1/(expo-1));      


U_new = tmp./(ones(cluster_n, 1)*sum(tmp));

mf = U_new.^expo;       
center = mf*data./((ones(size(data, 2), 1)*sum(mf'))');
