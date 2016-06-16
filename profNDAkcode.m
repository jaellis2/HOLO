function [ outphik ] = profNDAkcode(phik, data)

% data = struct('tau', tau, 'ir', IR, 'il', IL, 'ar', ang_r, 'aw', ang_w, ...
%      'xgrid', xgrid, 'sig_t', sig_t, 'sig_s', sig_s, 'sig_f', sig_f, 'q', q);

h = data.tau/(length(data.xgrid) - 1);


phi = phik(1:end-1);

if(h*sum(phi) ~= 1)
    phi = phi / (h * sum(phi));
end

k = phik(end);

Np_update = 1;

[phiHO, ephiHO, ejHO] = running_HOk(phi, k, data, Np_update, data.newfilename);

[phiHO, ephiHO, ejHO] = filterFlux(phiHO, ephiHO, ejHO, data);

[D] = d_calc3(ephiHO, ejHO, data);

phiHOb = [phiHO(1); phiHO(end)];

[ newphi ] = LOk_problem(ephiHO, k, phiHOb, D, data)


outphi = (newphi(1:end-1) + newphi(2:end))/2;

% outphiHO = (ephiHO(1:end-1) + ephiHO(2:end))/2;

normed_outphi = outphi / (h * sum(outphi));

outk = calc_k(outphi, phi, k, data);

% outphiHO = (ephiHO(1:end-1) + ephiHO(2:end))/2;

outphik = [normed_outphi; outk];

end

