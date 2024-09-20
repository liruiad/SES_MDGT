function x=DEE(c,L,N)
x=zeros(L,1);
  for k=0:L-1
      for n=0:N-1
          x(k+1)=x(k+1)+c(k+1,n+1)*exp(-j*2*pi*n*k/N);
      end
  end
x=real(x);
end

