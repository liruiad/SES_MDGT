function value = M_r(c,M,N)
value=0;
alpha=3;
for m=0:M-1;
    for n=0:N-1;
        value=value+c(m+1,n+1)^alpha;
    end
end
value=1/(1-alpha)*log2(value);
end
