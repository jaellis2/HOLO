function [ outpsi ] = Back_Sweep(phi, data)

%data = struct('tau', TAU, 'ir', IR, 'il', IL, 'w', W, 's', s, 'ar', ang_r, 'aw', ang_w, 'x', spat_grid);

na = length(data.ar);
nx = length(data.xgrid);
h = data.tau/(nx-1);

sig_t = data.st;
% sig_s = phi.*data.sig_s/2;

csig_s = (data.sig_s(1:end-1) + data.sig_s(2:end))/2;

% S = phi.*data.sig_s/2 + data.q_sour; %EDGE

S = phi.*(csig_s)/2 + (data.q)/2;

backward = [1:(na/2)]';


outpsi = zeros(nx,length(backward));

outpsi(nx,:) = data.ir*ones(1,length(backward));

st_avg = (sig_t(1:end-1) + sig_t(2:end))/2;


for i = (nx-1):-1:1
    
    p1 = 1./(-data.ar(backward) + (h*st_avg(i))/2);
%     p2 = (h*(sig_s(i+1)+sig_s(i))/2);
    p2 = h*S(i);

    p3 = (-data.ar(backward)-(h*st_avg(i))/2).*outpsi(i+1,:)';
    
    outpsi(i,:) = (p1.*(p2+p3))';
end



end
