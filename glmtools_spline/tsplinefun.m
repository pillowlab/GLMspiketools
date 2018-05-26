function [f,df,ddf] = tsplinefun(x,splfun,tfun)
%  [f,df,ddf] = tsplinefun(x,splfun,tfun)
%
%  Computes the function:
%     f(x) = tfun(splfun(x)) and its first and second derivatives
%  where tfun is a function handle


switch nargout
    case 0, 
	f = tfun(splfun(x));
    case 1,
	f = tfun(splfun(x));
    case 2,
	[g,dg] = splfun(x);
	[f,df] = tfun(g);
	df = df.*dg;
    case 3,
	[g,dg,ddg] = splfun(x);
	[f,dfg,ddfg] = tfun(g);
	df = dfg.*dg;
	ddf = ddfg.*dg.^2 + dfg.*ddg;
end
