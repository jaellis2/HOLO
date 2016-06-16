function [outphi, outephi, outJcurr] = filterFlux(phi, ephi, Jcurr, data)

regions = 6;
order = 3;

% span_v = 7;


xgrid = data.xgrid;

cgrid = (xgrid(1:end-1) + xgrid(2:end))/2;

tau = data.tau;
h = tau / (length(xgrid) - 1);

filter_x = [(tau/regions/2):(tau/regions):tau-(tau/regions/2)];




PP_phi = spap2(augknt([0, filter_x, tau], order), order, cgrid, phi);

PP_ephi = spap2(augknt([0, filter_x, tau], order), order, xgrid, ephi);

PP_Jcurr = spap2(augknt([0, filter_x, tau], order), order, xgrid, Jcurr);


outphi = fnval(PP_phi, cgrid);
outephi = fnval(PP_ephi, xgrid);
outJcurr = fnval(PP_Jcurr, xgrid);

% 
% outphi = smooth(phi, span_v);
% outephi = smooth(ephi, span_v);
% outJcurr  = smooth(Jcurr, span_v);


end