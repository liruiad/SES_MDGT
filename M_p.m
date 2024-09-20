function value = M_p(c,M,N)
value=0;
p=2;
for m=0:M-1;
    for n=0:N-1;
        value=value+c(m+1,n+1)^(1/p);
    end
end
value=value^p;
end
 
