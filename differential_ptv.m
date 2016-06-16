function [Pv] = differential_ptv(v, data)
%   phi's at centers

cn = length(v);
nx = length(data.xgrid);

h = data.tau/(nx-1);

sig_s = data.sig_s;
sig_t = data.st;

Pv = zeros(cn,1);


% if (isfield(data, 'xi'))
% xi = data.xi;
% phiHOv = xi*Fx_data.aw; % CELL EDGE: NX x 1
%     
%     cphiHOv = (phiHOv(1:end-1)+phiHOv(2:end))/2;
% 
% 
%     Pvb = v([1;end])-cphiHOv([1;end]);
% end
% 
% if(~isfield(data, 'xi'))
%     Pvb = [v(1); v(end)];
% 
% end

D = data.Dphi;

Pvb = [v(1); v(end)];


%%%%%%%%%FIND F(PHI)%%%%%%%%%%%



csig_t = (sig_t(1:end-1) + sig_t(2:end))/2;     %Total Cross-Sections at cell faces
csig_s = (sig_s(1:end-1) + sig_s(2:end))/2;     %Scattering Cross-Sections at cell faces




e = ones(cn-2,1);
Dmin = zeros(cn-2,1);
Dplus = zeros(cn-2,1);

et1 = (1./csig_t(3:(cn))).*e;
et2 = (1./csig_t(2:(cn-1))).*e;
et3 = (1./csig_t(1:(cn-2))).*e;

Dmin(1:(end-1)) = (1/(2*h))*D(2:(end-2));
Dplus(2:end) = (1/(2*h))*D(3:(end-1));

phimat = (-1/(3*h^2))*spdiags([et1 -2*et2 et3], -1:1, cn-2, cn-2);

phimat = phimat + spdiags([-Dmin.*e Dplus.*e], [-1,1], cn-2, cn-2); %

phimat = phimat + spdiags([(csig_t(2:(cn-1)) - csig_s(2:(cn-1))).*e] , 0, cn-2, cn-2);


b = v(2:end-1);

% b(1) = b(1) + (1/(3*csig_t(1)*h^2))*v(1);
% b(end) = b(end) + (1/(3*csig_t(end)*h^2))*v(end);


Pv(2:end-1) = phimat\b;
Pv(1) = Pvb(1);
Pv(end) = Pvb(2);

end