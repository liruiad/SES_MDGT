function c = DET(x,win,L,N)
 c=zeros(L,N);
 M0=L/N;
 for k=0:L-1;
   RR=zeros(N,1);
   for q=0:N-1;
      for p=0:M0-1;
          ii=mod(p*N+q,L);
          RR(q+1)=RR(q+1)+x(ii+1)*win(k+1,ii+1);
      end
    end
    c(k+1,:)=fft(RR)';
 end

end

