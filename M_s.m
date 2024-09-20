function value = M_s(c,M,N)
value=0;
for m=0:M-1;
    for n=0:N-1;
        if(abs(c(m+1,n+1))==0)
            c(m+1,n+1)=0.00000000001;
        end
        value=value+c(m+1,n+1)*log2(c(m+1,n+1));
    end
end 
value=-value;
end


