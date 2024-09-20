clc;
clear all;
L=256;
R=4;
sig=[5 10 16 30]; 
h=zeros(R,L);
for ix=0:R-1;
  for iy=0:L-1;
          sigma=sig(ix+1);
          h(ix+1,iy+1)=1/(pi*sigma^2)^0.5*exp(-(iy-0.5*(L-1))^2/(2*sigma^2));
  end
end

% save h ;
plot(h(1,:),'b');
hold on;
plot(h(2,:),'r');
hold on;
plot(h(3,:),'g');
hold on;
plot(h(4,:),'k');


N=128;
M0=L/N;
N0_p=[2 4 8 16];
M_p=L./N0_p;
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
figure;
plot(g0(1:L,1),'b');
hold on;
plot(g0((L+1):2*L),'r');
hold on;
plot(g0((2*L+1):3*L),'g');
hold on;
plot(g0((3*L+1):4*L),'k');

g=zeros(R,L);
for ix=0:R-1
    g(ix+1,:)=g0((ix*L+1):(ix+1)*L);
end
win=zeros(L,L);
for n=0:L-1;
    for l=0:L-1;
        for r=0:R-1;
            for m=0:M_p(r+1)-1;
                ix=mod(l-m*N0_p(r+1),L);
                iy=mod(n-m*N0_p(r+1),L);
                win(n+1,l+1)=win(n+1,l+1)+g(r+1,ix+1)*h(r+1,iy+1);
            end
        end
    end
end
figure_tvw(win(:,L/2));
