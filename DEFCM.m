function [v,u]=DEFCM(im1,im2,pre1,pre2,cluster_n,Lsw,pow,brain,par)

[row,col]=size(pre1);
ccen=[2,120,0.00001,0];
[v1,u1]=fcm(im1,cluster_n,ccen);
center=v1;
u_adapt=u1;
%%%%%%%%%%%%%%%%%%%%%%

ws=floor(Lsw/2);
epsilon=0.00001;
[im_row,im_col]=size(pre1);
w_orig=zeros(cluster_n,size(im1,1));
w_den=zeros(cluster_n,size(im1,1));
w_orig2=zeros(size(pre1));
w_den2=w_orig2;

for i=1:im_row
    for j=1:im_col
        c1=max(1,j-ws);
        c2=min(im_col,j+ws);
        r1=max(1,i-ws);
        r2=min(im_row,i+ws);
        Li2=pre1(r1:r2,c1:c2);
        Li_max=max(Li2(:));
        Li_den2=pre2(r1:r2,c1:c2);
        
        Li=reshape(Li2,[1,(c2-c1+1)*(r2-r1+1)]);
        
        sigma2Li=std(Li)^2;
        if  sigma2Li==0
            sigma2Li=epsilon;
        end
        w1_orig=exp((-1*((Li-pre1(i,j)).^2)/(sigma2Li)));
        
        w1_den=exp((-1*((Li-pre2(i,j)).^2)/(sigma2Li)));
        
        if  sum(w1_orig)+sum(w1_den)==0
            w_orig2(i,j)=0.5;
            w_den2(i,j)=0.5;
        else
            
            w_orig2(i,j)=(sum(w1_orig))/(sum(w1_orig)+sum(w1_den));
            w_den2(i,j)=(sum(w1_den))/(sum(w1_orig)+sum(w1_den));
        end
        
    end
end

hhh=conv2Dto1D(w_orig2);
hhh=hhh';

w_den(1,:)=(conv2Dto1D(w_den2))';

w_orig(1,:)=hhh;
for ppp=2:cluster_n
    w_orig(ppp,:)=w_orig(1,:);
    w_den(ppp,:)=w_den(1,:);
    
end

expo=2;

im1_n = size(im1, 1);
im2_n = size(im2, 1);
in_n = size(im1, 2);
max_iter=100;
display=1;
min_impro=0.00001;
cen2=[1000; 1000; 1000; 1000];
cen2=cen2(1:cluster_n);


upper_data=im1;
lower_data=im2;



center=center(1:cluster_n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Main loop
for i = 1:max_iter,
    
    [U1,ccccc,Kdist1] = mystepfcm(im1, cluster_n, expo,center);
    [U2,ccccc,Kdist2] = mystepfcm(im2, cluster_n, expo,center);
    
    
    if display,
        fprintf('Iteration count = %d\n', i);
    end
    
    upper=U2;
    lower=U1;
    
    mf1 = w_orig.*U1.^expo;
    mf1=mf1.*Kdist1;
    mf2 = w_den.* U2.^expo;
    mf2=mf2.*Kdist2;
    center = (mf1*im1+mf2*im2)./(((ones(size(im1, 2), 1)*sum(mf1'))')+((ones(size(im2, 2), 1)*sum(mf2'))'));
 
    if i > 1,
        if norm(cen2-center) < min_impro, break; end,
    end
    
    cen2=center;
    
    
end

u=(U1.*w_orig+U2.*w_den);
v=center;


