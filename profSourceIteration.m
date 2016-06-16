function [ outphiHO ] = profSourceIteration(phi, data)

% data = struct('tau', tau, 'ir', IR, 'il', IL, 'ar', ang_r, 'aw', ang_w, ...
%      'xgrid', xgrid, 'sig_t', sig_t, 'sig_s', sig_s, 'sig_f', sig_f, 'q', q);

% nx = length(data.xgrid);
% na = length(data.ar);
% 
% psi = zeros(nx,na);
% 
% psi(:,((na/2)+1):na) = Forward_Sweep(phi, data);
% psi(:,1:(na/2)) = Back_Sweep(phi, data);
% 
% eoutphi = psi*data.aw;
% 
% 
% coutphi = (eoutphi(1:end-1) + eoutphi(2:end))/2;

Np_update = 1;

[phiHO, ~, ~] = running_HO(phi, data, Np_update, data.newfilename);

% [phiHO, ~, ~] = filterFlux(phiHO, ephiHO, ejHO, data);

outphiHO = phiHO;

% outphiHO = (ephiHO(1:end-1) + ephiHO(2:end))/2;



end

