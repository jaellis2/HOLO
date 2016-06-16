function [ outphiHO ] = profNDASource(phi, data)

% data = struct('tau', tau, 'ir', IR, 'il', IL, 'ar', ang_r, 'aw', ang_w, ...
%      'xgrid', xgrid, 'sig_t', sig_t, 'sig_s', sig_s, 'sig_f', sig_f, 'q', q);

Np_update = 1;

[phiHO, ephiHO, ejHO] = running_HO(phi, data, Np_update, data.newfilename);


% [phiHO, ephiHO, ejHO] = filterFlux(phiHO, ephiHO, ejHO, data);


[D] = d_calc3(ephiHO, ejHO, data);

% phiHOb = [phiHO(1); phiHO(end)]; % ****

phiHOb = [ephiHO(1); ephiHO(end)];

[ newphi ] = LO_problem(phiHOb, D, data);


outphiHO = (newphi(1:end-1) + newphi(2:end))/2;

% outphiHO = (ephiHO(1:end-1) + ephiHO(2:end))/2;



end

