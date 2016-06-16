function [D,  phiHOb,  phiHO, jHO ] = d_calc(psi, data)
    
    nx = length(data.xgrid);
    h = data.tau/(nx-1);
    
    phiHO = psi*data.aw; % EDGES
        
    sig_t = data.st;            % Edges
 
    csig_t = data.csig_t;       % Centers
    
    
    phiHOb(1) = phiHO(1);
    phiHOb(2) = phiHO(nx);

    
    
    %%%%%%%%%%%%%%%%%
%     CELL FACES D/J %
    %%%%%%%%%%%%%%%%%
    
    mid_psi = (psi(2:nx,:)+psi(1:(nx-1),:))/2;     

    jHO = (mid_psi*(data.ar.*data.aw));     %Centers
    
         
%     D = (jHO + (1/(3*data.st))*(phiHO(2:nx) - phiHO(1:(nx-1)))/h);    % Good
    
%     D = (jHO + (1./(3./2*(sig_t(2:nx) + sig_t(1:(nx-1))))).*(phiHO(2:nx) - phiHO(1:(nx-1)))/h);    % Good
%     D = (2./(phiHO(2:nx) + phiHO(1:(nx-1)))).*D;
%     
    D = (jHO + (1./(3*(csig_t))).*(phiHO(2:nx) - phiHO(1:(nx-1)))/h);    % Good
    D = (2./(phiHO(2:nx) + phiHO(1:(nx-1)))).*D;
    

end

