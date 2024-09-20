clear all;clc ;
close all;
L=512;
fs=256;
x=[];
f=50;
for t=(1:1:L)/fs;
  if(t<1)
      x=[x 220*sqrt(2)*cos(2*pi*50*t)];
  end 
  if(t>=1)
      x=[x 220*sqrt(2)*(cos(2*pi*50*t)*(1+20/100*1/2*cos(2*pi*15*t)))];
  end  
end
x=mapminmax(x,-1,1);
x=x';
plot(x);
R=4;

sigma=[1 2 16 32];
for r=0:R-1;
    for l=0:L-1;
       sig=sigma(r+1);
        h(r+1,l+1)=1/(2*pi*sig^2)^0.5*exp(-(l-0.5*(L-1))^2/(2*sig^2)); 
    end
end

N=512;
M0=L/N;
N0_p=[1 4 4 1];
M=L./N0_p;
N0=max(N0_p);
g0=zeros(L*R,1);
H=[];
V=zeros(M0*N0,1);
V(1,1)=1/N;
for p=1:R
    H0=zeros(M0*N0,L);
   for m=0:M0-1
      for n=0:N0-1
         for k=0:L-1
           k1=mod(k+m*N,L);
           ix=mod(m*N0+n,M0*N0);
           if(n<N0_p(p))
             H0(ix+1,k+1)=1/N0_p(p)*h(p,k1+1)*exp(j*2*pi*n*k/N0_p(p));  
           else
             H0(ix+1,k+1)=1/N0_p(p)*h(p,k1+1)*0;  
           end
         end
      end
   end
   H=[H H0];
   clear H0;
end

G=inv(double(H*H'));
g0=real((H'*G)*V);
g=zeros(R,L);
for ix=0:R-1
    g(ix+1,:)=g0((ix*L+1):(ix+1)*L);
end
% figure;
% plot(gamma(1,:));hold on;
% plot(gamma(2,:));hold on;
% plot(gamma(3,:));hold on;
% plot(gamma(4,:));hold on;
% % title('analysis windows');
win=zeros(L,L);
for n=0:L-1;
    for l=0:L-1;
        for r=0:R-1;
            for m=0:M(r+1)-1;
                ix=mod(l-m*N0_p(r+1),L);
                iy=mod(n-m*N0_p(r+1),L);
                win(n+1,l+1)=win(n+1,l+1)+g(r+1,ix+1)*h(r+1,iy+1);
            end
        end
    end
end
figure; 
plot(win(:,256));

c=DET(x,win,L,N);
E0=sum(sum(abs(c).^2));
cv0=abs(c).^2/E0;
value_p = M_p(cv0,L,N);
value_s = M_s(cv0,L,N);
value_r = M_r(cv0,L,N);
value_nr= M_nr(cv0,L,N);
v_es=[value_p value_s value_r value_nr]
createfigureEX3(abs(c(:,1:N/2)').^2);
iter=10;
mu=0.5;
gamma=0.8;
ct=c;
K=floor(L*N*0.02);
it=0;
while(it<iter)
    xt=DEE(ct,L,N);
    xt_hat=x-xt;
    ct_hat=DET(xt_hat,win,L,N);
    ct_bar=ct+mu*ct_hat;
    c_bar=reshape(abs(ct_bar),L*N,1);
    c_bar=sort(c_bar,'descend');
    lambda=sqrt(96)/(9*mu*(1+mu*gamma)^0.5)*c_bar(K+1)^1.5;
    ct=shrinkge(lambda,mu,gamma,ct_bar);
    it=it+1;
end
createfigureEX3(abs(ct(:,1:N/2)').^2);
E1=sum(sum(abs(ct).^2));
cv1=abs(ct).^2/E1;
value_p = M_p(cv1,L,N);
value_s = M_s(cv1,L,N);
value_r = M_r(cv1,L,N);
value_nr= M_nr(cv1,L,N);
v_ses=[value_p value_s value_r value_nr]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M1=L/N0_p(1);
c=zeros(M1,N);
c=SDGT(L,M1,N,x,h(1,:)');
tt=c(1:M1/2,:);
c(1:M1/2,:)=c(M1/2+1:M1,:);
c(M1/2+1:M1,:)=tt;
createfigureEX3(abs(c(:,1:N/2)').^2);
E2=sum(sum(abs(c).^2));
cv2=abs(c).^2/E2;
value_p = M_p(cv2,M1,N);
value_s = M_s(cv2,M1,N);
value_r = M_r(cv2,M1,N);
value_nr= M_nr(cv2,M1,N);
v_1=[value_p value_s value_r value_nr]

M1=L/N0_p(2);
c=SDGT(L,M1,N,x,h(2,:)');
tt=c(1:M1/2,:);
c(1:M1/2,:)=c(M1/2+1:M1,:);
c(M1/2+1:M1,:)=tt;
createfigureEX3(abs(c(:,1:N/2)').^2);
E2=sum(sum(abs(c).^2));
cv2=abs(c).^2/E2;
value_p = M_p(cv2,M1,N);
value_s = M_s(cv2,M1,N);
value_r = M_r(cv2,M1,N);
value_nr= M_nr(cv2,M1,N);
v_2=[value_p value_s value_r value_nr]

M1=L/N0_p(3);
c=SDGT(L,M1,N,x,h(3,:)');
tt=c(1:M1/2,:);
c(1:M1/2,:)=c(M1/2+1:M1,:);
c(M1/2+1:M1,:)=tt;
createfigureEX3(abs(c(:,1:N/2)').^2);
E2=sum(sum(abs(c).^2));
cv2=abs(c).^2/E2;
value_p = M_p(cv2,M1,N);
value_s = M_s(cv2,M1,N);
value_r = M_r(cv2,M1,N);
value_nr= M_nr(cv2,M1,N);
v_3=[value_p value_s value_r value_nr]

M1=L/N0_p(4);
c=SDGT(L,M1,N,x,h(4,:)');
tt=c(1:M1/2,:);
c(1:M1/2,:)=c(M1/2+1:M1,:);
c(M1/2+1:M1,:)=tt;
createfigureEX3(abs(c(:,1:N/2)').^2);
E2=sum(sum(abs(c).^2));
cv2=abs(c).^2/E2;
value_p = M_p(cv2,M1,N);
value_s = M_s(cv2,M1,N);
value_r = M_r(cv2,M1,N);
value_nr= M_nr(cv2,M1,N);
v_4=[value_p value_s value_r value_nr]

M1=L/N0_p(5);
c=SDGT(L,M1,N,x,h(5,:)');
tt=c(1:M1/2,:);
c(1:M1/2,:)=c(M1/2+1:M1,:);
c(M1/2+1:M1,:)=tt;
createfigureEX3(abs(c(:,1:N/2)').^2);
E2=sum(sum(abs(c).^2));
cv2=abs(c).^2/E2;
value_p = M_p(cv2,M1,N);
value_s = M_s(cv2,M1,N);
value_r = M_r(cv2,M1,N);
value_nr= M_nr(cv2,M1,N);
v_5=[value_p value_s value_r value_nr]


