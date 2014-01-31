function rtn = ParLoP( S,N )
%PARLOP Summary of this function goes here
%   Detailed explanation goes here

rtn(1)=-88;
parfor i=1:length(S);
	
    try
	if S(i) > log2(N(i)),
		rtn(i) = -99;
	else
        	rtn(i) = getLimit(S(i),N(i));
	end
    catch exception
	rtn(i) = -88;
    end
end
end

