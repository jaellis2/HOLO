function outk = calc_k(outphi, phi, k, data)

sig_f = data.sig_f;
nu_f = data.nu_f;

h = data.tau/(length(data.xgrid) - 1);

outk = ((h * sum(sig_f .* nu_f .* outphi)) ./ (h * sum(sig_f .* nu_f .* phi))) * k;

end