function [ Jv ] = transport_jv(v, Fx_data)


%%%%%%
%   Profugus VERSION
%


nx = length(Fx_data.xgrid);
na = length(Fx_data.ar);
h = Fx_data.tau/(nx-1);
sig_t = Fx_data.st;
sig_s = Fx_data.sig_s;
phi = Fx_data.phiLO;   %CENTER: from the knl construction of data =  ['fx', fx, 'x', x, 'static_data', static_data]?

Fphi = Fx_data.Fphi;

Dphi = Fx_data.Dphi;
phiHO = Fx_data.phiHO;
jHO = Fx_data.jHO;

cn = length(phi);
inn = cn-2;

csig_t = (sig_t(1:end-1) + sig_t(2:end))/2;
csig_s = (sig_s(1:end-1) + sig_s(2:end))/2;



Jv = zeros(nx-1,1); % for the Phi's at cell-centers

Lv = zeros(nx-1,1);
Nv = zeros(nx-1,1);

%   e = ones(nx-2,1);

%%%% D Operator %%%%
e1 = ones(inn+2,1);
e2 = ones(inn,1);

cddxc = spdiags([-e2 e2], [-1 1], inn, inn);
cd2dx2c = spdiags([e2 -2*e2 e2], [-1:1], inn, inn);

eddxc = spdiags([-e1 e1], [0 1], nx-1, nx);


cddxc = 1/(2*h).*cddxc;  % CENTER TO CENTER D/DX OPERATOR
cd2dx2c = 1/(h^2).*cd2dx2c; % CENTER TO CENTER D^2/DX^2 OPERATOR

eddxc = 1/(h)*eddxc; % EDGE TO CENTER D/DX OPERATOR


%%%%

%
%     dvdx = zeros(nx-1,1);
%     d2vdx2 = zeros(nx-1,1);
%
%      dvdx = cddxc*v;
%      dvdx(1) = dvdx(1) - v(1)/(2*h);
%      dvdx(end) = dvdx(end) + v(end)/(2*h);

%     d2vdx2 = cd2dx2c*v;
%     d2vdx2(1) = d2vdx2(1) + v(1)/(h^2);
%     d2vdx2(end) = d2vdx2(end) + v(end)/(h^2);
%
%      dsig_t = eddxc*(-1./(3*sig_t));
%     dsig_t(1) = dsig_t(1) + 1./(3*sig_t(1));
%     dsig_t(end) = dsig_t(end) - 1./(3*sig_t(end));


%     Lv = (-1./(3*sig_t)).*(d2vdx2) + (dvdx.*dsig_t) + (sig_t - sig_s).*v;      %Should be fine.



%Lv = (-1./(3*csig_t)).*(d2vdx2) + (dvdx.*dsig_t) + (csig_t - csig_s).*v;

%  Lv = (-1./(3*csig_t(2:end-1))).*(d2vdx2) + (csig_t(2:end-1) - csig_s(2:end-1)).*v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(inn,1);


%Questionable!

%     et1 = (1./csig_t(3:(inn))).*e;
%     et2 = (1./csig_t(2:(inn-1))).*e;
%     et3 = (1./csig_t(1:(inn-2))).*e;

et_avg = (1./csig_t(2:(end-1))).*e;


phimat = (-1/(3*h^2)).*spdiags([1*(et_avg) -2*(et_avg) 1*(et_avg)], -1:1, inn, inn);

phimat = phimat + spdiags([(csig_t(2:(end-1)) - csig_s(2:(end-1))).*e] , 0, inn, inn);

% b(1) = b(1) + (1/(3*csig_t(1)*h^2))*v(1);
% b(end) = b(end) + (1/(3*csig_t(end)*h^2))*v(end);


Lv = phimat*v(2:end-1);

Lv(1) = Lv(1) + (-1/(3*csig_t(2)*h^2))*v(1);
Lv(end) = Lv(end) + (-1/(3*csig_t(end-1)*h^2))*v(end);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Nonlinear Operator%%%%

% xi = zeros(nx,na);

% HO psi

vs = [v];

newdata = Fx_data;


newdata.q = zeros(length(v),1);
newdata.ir = 0;
newdata.il = 0;

% newdata.sig_s = zeros(inn+3,1);



% xi(:,((na/2)+1):na) = Forward_Sweep(vs, newdata);    % Good
% xi(:,1:(na/2)) = Back_Sweep(vs, newdata);            % Good


%%%%%%
%
%

Np_update = 1; % Scales number of particles down to equal total #


%%%%%
%   Split Neg/Pos
%

source_p = vs;
source_p(source_p <= 0) = 0;

source_n = -vs;
source_n(source_n <= 0 ) = 0;

norm_p = h * norm(source_p, 1);

norm_n = h * norm(source_n, 1);


cphiHOv_p = 0;
phiHOv_p = 0;
jHOv_p = 0;
cphiHOv_n = 0;
phiHOv_n = 0;
jHOv_n = 0;

if(norm_p == 0)
    
    newdata.np_factor = 1;
    
    [cphiHOv_n, phiHOv_n, jHOv_n] = running_HO(source_n, newdata, Np_update, newdata.newfilename); % CENTERS, EDGES
    
    [cphiHOv_n, phiHOv_n, jHOv_n] = filterFlux(cphiHOv_n, phiHOv_n, jHOv_n, newdata);

    
    
elseif(norm_n == 0)
    
    newdata.np_factor = 1;
    
    [cphiHOv_p, phiHOv_p, jHOv_p] = running_HO(source_p, newdata, Np_update, newdata.newfilename); % CENTERS, EDGES
    
    
else
    
    np_factor_p = round(norm_p / (norm_p + norm_n), 4);
    np_factor_n = round((1 / np_factor_p) * norm_n / (norm_p + norm_n), 4);
    
    final_np_factor = round(1 / (np_factor_p) * 1 / (np_factor_n), 4);
    
    
    
    newdata.np_factor = np_factor_p;
    
    [cphiHOv_p, phiHOv_p, jHOv_p] = running_HO(source_p, newdata, Np_update, newdata.newfilename); % CENTERS, EDGES
    
    [cphiHOv_p, phiHOv_p, jHOv_p] = filterFlux(cphiHOv_p, phiHOv_p, jHOv_p, newdata);
    
    

    newdata.np_factor = np_factor_n;
    
    [cphiHOv_n, phiHOv_n, jHOv_n] = running_HO(source_n, newdata, Np_update, newdata.newfilename); % CENTERS, EDGES    
    
    [cphiHOv_n, phiHOv_n, jHOv_n] = filterFlux(cphiHOv_n, phiHOv_n, jHOv_n, newdata);

    
    
    newdata.np_factor = final_np_factor;
    
    worked = addSource(newdata.newfilename, vs, Np_update, [1; Fx_data.Np], newdata);
    
end


cphiHOv = cphiHOv_p - cphiHOv_n;
phiHOv = phiHOv_p - phiHOv_n;
jHOv = jHOv_p - jHOv_n;






%%%%%
% Negative Profugus
%


% [cphiHOv, phiHOv, jHOv] = running_HO(vs, newdata, Np_update, newdata.newfilename); % CENTERS, EDGES




%%%%%%%%%%%%%%%%%%%%%%%%


% phiHOv = xi*Fx_data.aw; % CELL EDGE: NX x 1
% jHOv = (xi*(Fx_data.ar.*Fx_data.aw)); % CELL EDGE: NX x 1
%NOT INTERIOR NOT CENTER

% cphiHOv = (phiHOv(1:end-1) + phiHOv(2:end))/2; % CELL FACE    % 1st ELEMENT
cjHOv = (jHOv(1:end-1) + jHOv(2:end))/2;
%NOT INTERIOR






cphiHO = (phiHO(1:end-1) + phiHO(2:end))/2;
% cjHO = (jHO(1:end-1) + jHO(2:end))/2;
cjHO = (jHO(2:end) + jHO(1:end-1))/2;

dphiHOv = eddxc*phiHOv; %EXACT/GOOD %FACES            %2nd ELEMENT -- SOLVED



dphiHO = eddxc*phiHO; %EXACT/GOOD   %FACES
%NOT INTERIOR

dDv = cphiHO.*(cjHOv + 1./(3*csig_t).*dphiHOv);
dDv = dDv - (cjHO + 1./(3*csig_t).*dphiHO).*cphiHOv;
dDv = dDv./((cphiHO).^2);
%NOT INTERIOR

% Dphi = (cjHO + 1./(3*csig_t).*dphiHO)./cphiHO;
%NOT INTERIOR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Jeff's Thesis ~
%
% dDv1 = cphiHO.*(cjHOv + 1./(3*csig_t).*dphiHOv)./(cphiHO.^2);
% dDv2 = (cphiHOv.*Dphi)./(cphiHO);
%
% dDv_new = dDv1 - dDv2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%INTERIOR%%%%
%     Nv = cddxc*(dDv(2:end-1).*phi(2:end-1) + Dphi(2:end-1).*v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Product Rule

% Nv_a = (cddxc*(dDv(2:end-1))).*phi(2:end-1);
% Nv_a(1) = Nv_a(1) - dDv(1)/(2*h)*phi(2);
% Nv_a(end) = Nv_a(end) + dDv(end)/(2*h)*phi(end-1);
%
% Nv_b = dDv(2:end-1).*(cddxc*(phi(2:end-1)));
% Nv_b(1) = Nv_b(1) - dDv(2)*phi(1)/(2*h);
% Nv_b(end) = Nv_b(end) + dDv(end-1)*phi(end)/(2*h);
%
% Nv_c = (cddxc*(Dphi(2:end-1))).*v(2:end-1);
% Nv_c(1) = Nv_c(1) - dDv(1)/(2*h)*v(2);
% Nv_c(end) = Nv_c(end) + dDv(end)/(2*h)*v(end-1);
%
% Nv_d = Dphi(2:end-1).*(cddxc*(v(2:end-1)));
% Nv_d(1) = Nv_d(1) - Dphi(2)*v(1)/(2*h);
% Nv_d(end) = Nv_d(end) + Dphi(end-1)*v(end)/(2*h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%USE EXTERIOR BOUNDARIES%%%%
%     Nv2 = cddxc*(dDv(2:end-1).*phi(2:end-1)) + cddxc*(Dphi(2:end-1).*v(2:end-1));
%     nva = (1/(2*h))*(dDv(1)*phi(1) + Dphi(1)*v(1));
%     nvb = (1/(2*h))*(dDv(end)*phi(end) + Dphi(end)*v(end));
%
%     Nv2(1) = Nv2(1) - nva;
%     Nv2(end) = Nv2(end) + nvb;

% %

% Nv = Nv_a + Nv_b + Nv_c + Nv_d;


Nv3_a = cddxc*(dDv(2:end-1).*phi(2:end-1));
Nv3_a(1) = Nv3_a(1) - (1/(2*h))*(dDv(1)*phi(1));
Nv3_a(end) = Nv3_a(end) + (1/(2*h))*(dDv(end)*phi(end));

Nv3_b = cddxc*(Dphi(2:end-1).*v(2:end-1));
Nv3_b(1) = Nv3_b(1) - (1/(2*h))*(Dphi(1)*v(1));
Nv3_b(end) = Nv3_b(end) + (1/(2*h))*(Dphi(end)*v(end));

Nv3 = Nv3_a + Nv3_b;

Jvb = zeros(cn,1);

Jvb(:) = v - cphiHOv;


%%%%%%%%%%%


Jv(2:end-1) = (Lv + Nv3); % INTERIOR CELL FACES

Jv(1) = Jvb(1);
Jv(end) = Jvb(end);

end