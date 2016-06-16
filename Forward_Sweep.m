function [ outpsi ] = Forward_Sweep(phi, data)

%data = struct('tau', TAU, 'ir', IR, 'il', IL, 'w', W, 's', s, 'ar', ang_r, 'aw', ang_w, 'x', spat_grid);


na = length(data.ar);
nx = length(data.xgrid);
h = data.tau/(nx-1);


forward = [(na/2+1):na]';


outpsi = zeros(nx,length(forward));

sig_t = data.st;

csig_s = (data.sig_s(1:end-1) + data.sig_s(2:end))/2;


% S = phi.*data.sig_s/2 + data.q_sour; %EDGE

S = phi.*(csig_s)/2 + (data.q)/2;

outpsi(1,:) = data.il*ones(1,length(forward)); % check

st_avg = (sig_t(1:end-1) + sig_t(2:end))/2;


for i = 2:(nx)

    
    p1 = (1./(data.ar(forward) + (h*st_avg(i-1))/2));
%     p2 = (h*(S(i)+S(i-1))/2);
    p2 = h*S(i-1);
    p3 = (data.ar(forward)-(h*st_avg(i-1))/2).*outpsi(i-1,:)';
    
    outpsi(i,:) = (p1.*(p2+p3))';
end




end

