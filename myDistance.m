function [D1]=myDistance(V1,c,I1)

[m,n]=size(I1);

D1=zeros(c,m,n);

I_mean=mean(I1(:));

ff=(I1-I_mean).^2;
flag=sum(ff(:));

hh=((I1-I_mean).^2-flag/(m*n)).^2;
h=sum(hh(:));

h=1/sqrt(h/(m*n-1));
for k=1:c
    D1(k,:,:)=exp((-1)*h*(I1-V1(k)).^2);
end

            
            
            