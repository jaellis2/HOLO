function [ newphi ] = LO_problem(phiHOb, D, data)


    nx = length(data.xgrid);
    h = data.tau/(nx-1);
    
    sig_s = data.sig_s;
    sig_t = data.st;
    
    csig_t = data.csig_t;
    csig_s = data.csig_s;


    newphi = zeros(nx,1);
    


    %%%%%%%%%%%%%%%%%%%
%     CELL EDGE MATRIX %
    %%%%%%%%%%%%%%%%%%%
%     
%     e = ones(nx-2,1);
%     Dmid1 = zeros(nx-2,1);
%     Dmid2 = zeros(nx-2,1);
% 
%     
%     phimat = (-1/(3*data.st*h^2))*spdiags([e -2*e e], -1:1, nx-2, nx-2);
%     
%     
%     Dmid1(1:end-1) = (1/(2*h))*(D(2:(nx-2)));
%     Dmid2(2:end) = (1/(2*h))*(D(3:(nx-1)));
%     phimat = phimat + spdiags([-Dmid1.*e Dmid2.*e], [-1 1], nx-2, nx-2);
% 
%     
%     
%     
%     phimat = phimat + spdiags([(data.st - sig_s(2:(nx-1))).*e], 0, nx-2, nx-2);
%     
%    
%     b = (data.q_sour*ones(nx-2,1));
%     b(1) = b(1) + ((1/(3*data.st*h^2)) + (D(1)/(2*h)))*phiHOb(1);
%     b(end) = b(end) + ((1/(3*data.st*h^2)) - (D(nx)/(2*h)))*phiHOb(2);
%     
%     
%     newphi(2:(nx-1)) = phimat\b;
%     newphi(1) = phiHOb(1);
%     newphi(nx) = phiHOb(2);

    
    %%%%%%%%%%%%%%%%%%%
%     CELL FACE MATRIX %
    %%%%%%%%%%%%%%%%%%%
    
    e = ones(nx-2,1);
    Dmin = zeros(nx-2,1);
    Dmid = zeros(nx-2,1);
    Dplus = zeros(nx-2,1);
    etmin = zeros(nx-2,1);
    etplus = zeros(nx-2,1);
    
%     et1 = (1./sig_t(3:(nx))).*e;
%     et2 = (1./sig_t(2:(nx-1))).*e;
%     et3 = (1./sig_t(1:(nx-2))).*e;
% 
%     phimat = (-1/(3*h^2))*spdiags([et1 -2*et2 et3], -1:1, nx-2, nx-2);
%     
    
    
    etmin(1:end) = 1./(csig_t(1:end-1));
    etplus(1:end) = 1./(csig_t(2:end));


    phimat = (-1/(3*h^2))*spdiags([etmin -(etmin+etplus) etplus], -1:1, nx-2, nx-2);


   
    
    Dmin(1:(end-1)) = (1/(2*h))*D(2:(nx-2));
    
    Dmid(1:end) = (1/(2*h))*(D(2:(nx-1)) - D(1:(nx-2)));
    
    Dplus(2:end) = (1/(2*h))*D(2:(nx-2));

    phimat = phimat + spdiags([-Dmin.*e Dmid.*e Dplus.*e], -1:1, nx-2, nx-2);
    
    phimat = phimat + spdiags([(sig_t(2:(nx-1)) - sig_s(2:(nx-1))).*e] , 0, nx-2, nx-2);
    
    
    
%     b = (data.q(2:end-1)); % ***

    b = data.q_sour*ones(nx-2, 1);
    
    b(1) = b(1) + ((1/(3*csig_t(1)*h^2)) + (D(1)/(2*h)))*phiHOb(1);
    b(end) = b(end) + ((1/(3*csig_t(end)*h^2)) - (D(nx-1)/(2*h)))*phiHOb(2);
    
    
    newphi(2:(nx-1)) = phimat\b;
    newphi(1) = phiHOb(1);
    newphi(nx) = phiHOb(2);


end
