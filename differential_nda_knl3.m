function [Fphi, Fx_data] = differential_nda_knl3(phi, data)
%   phi's at centers

%     data = struct('tau', TAU, 'ir', IR, 'il', IL, 'w', W, 's', S, ...
%             'ar', ang_r, 'aw', ang_w, 'xgrid', spat_grid, 'st', sig_t, 'sig_s', sig_s, 'q_sour', qsource);

data;


cn = length(phi);
nx = length(data.xgrid);
% na = length(data.ar);

h = data.tau/(nx-1);

sig_s = data.sig_s; %EDGES
sig_t = data.st;
% 
% psi = zeros(nx,na);
 Fphi = zeros(cn,1); %CENTERS
% 
% % HO psi
% psi(:,((na/2)+1):na) = Forward_Sweep(phi, data);    % Good
% psi(:,1:(na/2)) = Back_Sweep(phi, data);            % Good EDGES


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HO PROBLEM:  RUNNING PROFUGUS 
%



Np_update = 1;

[phiHO_raw, ephiHO_raw, jHO_raw] = running_HO(phi, data, Np_update, data.newfilename); % CENTERS, EDGES


[phiHO, ephiHO, jHO] = filterFlux(phiHO_raw, ephiHO_raw, jHO_raw, data);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate fudge factor
[D] = d_calc3(ephiHO, jHO, data); 

%%%%%%%%%FIND F(PHI)%%%%%%%%%%%


%Questionable
% Fb = [phi(1)-(phiHO(1)+phiHO(2))/2; phi(end)-(phiHO(end)+phiHO(end-1))/2];

Fb = [phi(1)-phiHO(1); phi(end)-phiHO(end)];


csig_t = (sig_t(1:end-1) + sig_t(2:end))/2;     %Total Cross-Sections at cell faces
csig_s = (sig_s(1:end-1) + sig_s(2:end))/2;     %Scattering Cross-Sections at cell faces



e = ones(cn-2,1);
Dmin = zeros(cn-2,1);
Dplus = zeros(cn-2,1);




et1 = (1./csig_t(3:(cn))).*e;
et2 = (1./csig_t(2:(cn-1))).*e;
et3 = (1./csig_t(1:(cn-2))).*e;

phimat = (-1/(3*h^2))*spdiags([et1 -2*et2 et3], -1:1, cn-2, cn-2);






Dmin(1:(end-1)) = (1/(2*h))*D(2:(cn-2));

Dplus(2:end) = (1/(2*h))*D(3:(cn-1));


phimat = phimat + spdiags([-Dmin.*e Dplus.*e], [-1,1], cn-2, cn-2);

phimat = phimat + spdiags([(csig_t(2:(cn-1)) - csig_s(2:(cn-1))).*e] , 0, cn-2, cn-2);



b = (data.q_sour*ones(cn-2,1)); % FIXED SOURCE, NO FISSION

b(1) = b(1) + ((1/(3*csig_t(1)*h^2)) + (D(1)/(2*h)))*phi(1);
b(end) = b(end) + ((1/(3*csig_t(end)*h^2)) - (D(end)/(2*h)))*phi(end);


%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%

Fphi(2:end-1) = phimat*phi(2:end-1);
Fphi(2:end-1) = Fphi(2:end-1) - b;
% Fphi(2) = Fphi(2) - (1/(3*csig_t(1)*h^2))*phi(1) - (D(1)/(2*h))*phi(1);
% Fphi(end-1) = Fphi(end-1) - (1/(3*csig_t(end)*h^2))*phi(end) + (D(end)/(2*h))*phi(end);
Fphi(1) = Fb(1);
Fphi(end) = Fb(2);


Fx_data = data;

Fx_data.Dphi = D;
Fx_data.phiHO = ephiHO;
Fx_data.jHO = jHO;
Fx_data.Fphi = Fphi;
Fx_data.phiLO = phi;

end