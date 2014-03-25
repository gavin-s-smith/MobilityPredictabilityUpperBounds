function rtn = getLimit( S,N )

syms x real;

E = (-(x*log2(x)+(1-x)*log2(1-x))+(1-x)*log2(N-1) ) - S;

rtn = vpasolve(E,0.9);

end

