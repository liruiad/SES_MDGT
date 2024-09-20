function value = M_nr(c,M,N)
v1=0;
v2=0;
alpha=0.8;%0.9;
beta=1.2; 
for m=0:M-1;
    for n=0:N-1;     
     v1=v1+c(m+1,n+1)^beta;
     v2=v2+c(m+1,n+1)^alpha;
    end
end
value=v1^(1/beta)/v2^(1/alpha);
end 
