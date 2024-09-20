clc;clear all;
L=512;N=L;
load('tpd_ex3.mat')
x=x;

TF=[];
% K=round(sqrt(1/alpha));
h=tftb_window(257,'Hanning',0.005);
% tic
% for k=0:L/N-1;
%     pos=1+k*N:N+k*N;
%     y=x(pos);
%     [tf,t,f]=tfrstft(y,1:N,N,h);
%     TF=[TF tf];
% end
% toc
[TF,t,f]=tfrstft(x,1:N,N,h);
createfigureEX3(abs(TF(1:N/2,:)).^2);
c1=TF';
E=sum(sum(abs(c1).^2));
S=abs(c1).^2/E;
M=L;
value_p=M_p(abs(S),L,N);
value_s=M_s(abs(S),L,N);
value_r=M_r(abs(S),L,N);
value_nr=M_nr(abs(S),L,N);
v_stft=[value_p value_s value_r value_nr]