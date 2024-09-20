clc;clear all;
L=512;N=L;
load('tpd_ex3.mat')
x=x;
% tic
% for k=0:Ls/N-1;
%     pos=1+k*N:N+k*N;
%     y=x(pos);
%     [s,t,f,E]=st(y,fs);
%     TF=[TF s];
% end
% toc
TF=st(x,256);
createfigureEX3(abs(TF(1:N/2,:)).^2);
c1=TF';
E=sum(sum(abs(c1).^2));
S=abs(c1).^2/E;
M=L;
value_p=M_p(abs(S),L,N/2);
value_s=M_s(abs(S),L,N/2);
value_r=M_r(abs(S),L,N/2);
value_nr=M_nr(abs(S),L,N/2);
v_st=[value_p value_s value_r value_nr]