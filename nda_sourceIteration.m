function [ coutphi, phiHO, jHO, D ] = nda_sourceIteration(phi, data)


nx = length(data.xgrid);
na = length(data.ar);


psi = zeros(nx,na);

% HO psi 
psi(:,((na/2)+1):na) = Forward_Sweep(phi, data);    % Good
psi(:,1:(na/2)) = Back_Sweep(phi, data);            % Good

% Calculate fudge factor
[D, phiHOb, phiHO, jHO] = d_calc(psi, data); %Good

eoutphi = LO_problem(phiHOb, D, data); 
coutphi = (eoutphi(1:end-1) + eoutphi(2:end))/2;


% 
% [D, phiHOb, phiHO, jHO] = d_calc2(psi, data); %Good
% 
% % Solve LO 
% coutphi = LO_problem2(phiHOb, D, data); 
% 

end

