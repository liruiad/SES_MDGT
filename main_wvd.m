clc;clear all;
L=512;N=256;
load('tpd_ex3.mat')
x=x;
x=hilbert(x);
[SWD,t,f] = tfrwv(x,1:L,N);
% SWD=zeros(size(WD));
%  nk=(L-64)/64+1;
%  for ix=0:nk-1;
%      for iy=0:nk-1;
%          pos1=(1:64)+ix*64;
%          pos2=(1:64)+iy*64;
%          wd=WD(pos1,pos2);
%          wd1= CS_WD(wd);
%          wd1=reshape(wd1,64,64);
%         SWD(pos1,pos2)=(SWD(pos1,pos2)+wd1);
%      end
%  end
% cc=SWD;
% min0=min(min(abs(cc)));
% max0=max(max(abs(cc)));
createfigureEX3(abs(SWD(1:N,:)).^2);
c1=SWD';
E=sum(sum(abs(c1).^2));
cv=abs(c1).^2/E;
% % value_s= M_s(cv,N,M)
value_p = M_p(cv,L,N);
value_s = M_s(cv,L,N);
value_r = M_r(cv,L,N);
value_nr= M_nr(cv,L,N);
v_wvd=[value_p value_s value_r value_nr]
