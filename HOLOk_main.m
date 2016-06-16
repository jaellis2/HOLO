%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HO/LO Profugus & NDA connectivity.
%
%   J. Austin Ellis
%   18 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

%%%%%%%%%
%
%   CHANGE TO PERSONAL FILE LOCATIONS
%
%


input_filename = 'my_mcfixedsrc_cart50.xml';
running_filename = 'running_mcfixedsrc_cart.xml';

xmc_pwd = '~/Documents/SOFTWARE/ProfugusB/ProfRelease/bin/';
HO_sol_pwd = '~/Documents/SOFTWARE/HOLO/HighOrderSolutions/';


np_factor = 1;

cs_factor = .2;
cf_factor = .5;

nu_f = 1.2;
chi = 1;

q_external = 0;

% c_factor = .8;

% c_factor = .9;
% 
% sig_f = .5;
% nu_f = 1.2;
% chi = 1
% sig_s = .2
% sig_t = 1





% holo_elem = {'x_edges'; 'general_source'; 'xs library'};
holo_elem = {'x_edges'; 'xs library'; 'Np'; 'problem_name'};
% mat_elem = {'sigma_t'; 'sigma_s0', 'nu_sigma_f', 'sigma_f', 'chi'};
mat_elem = {'sigma_t'; 'sigma_s0'};


input_data = parseXML(input_filename, running_filename, holo_elem, mat_elem);

%%%%%%%%
%
%  Set up parameters
%

tau = input_data.tau;
sig_t = input_data.sig_t;
sig_s = sig_t * cs_factor;
sig_f = sig_t * cf_factor;
Np = input_data.Np;
prob_name = input_data.prob_name;


xgrid = input_data.xgrid;
cgrid = (xgrid(1:end-1) + xgrid(2:end))/2.0;

NX = length(xgrid);


%%%%
%   HARD CODED INFO
%

sig_t = sig_t*ones(NX,1);
sig_s = sig_s*ones(NX,1);
csig_t = (sig_t(1:end-1) + sig_t(2:end))/2;
csig_s = (sig_s(1:end-1) + sig_s(2:end))/2;


qvec = q_external*ones(NX-1,1);

% qvec = zeros(NX-1,1);

IR = 0;
IL = 0;
W = 0;
S = 0;
angles = 20;

[ang_r, ang_w] = gauss(angles);

%%%%%

data = struct('tau', tau, 'ir', IR, 'il', IL, 'w', W, 's', S, 'Np', Np, ...
    'ar', ang_r, 'aw', ang_w, 'xgrid', xgrid, 'st', sig_t, ...
    'sig_s', sig_s, 'cgrid', cgrid, 'csig_t', csig_t, 'csig_s', csig_s, ...
    'sig_f', sig_f, 'nu_f', nu_f, 'chi', chi, ...
    'q', qvec, 'q_sour', q_external, 'prob_name', prob_name, 'filename', input_filename, ...
    'newfilename', running_filename, 'xmc_pwd', xmc_pwd, 'HO_sol_pwd', HO_sol_pwd, ...
    'np_factor', np_factor);

%
%  End Set-Up
%
%%%%%%%%%%

runProfugusSource = 0;
runProfugusAnderson = 1;
runProfugusNDA = 1;
runProfugusNDAAnderson = 1;
runProfugusNDA_KNL = 0;
runAndersonSource = 0;
runAndersonNDA = 0;

mlevel = 2;
maxit = 10;
atol = 1.d-5;
rtol = 1.d-5;

%Np0 = 20K or 200K

% 5mil particles  = 16sec. approx


%%%%%%%%%%
%
%  Begin Running
%

sourcef0 = ones(NX-1,1);     % CELL-CENTERED SCALAR FLUX INITIAL ITERATE

k0 = sig_f*nu_f;

source_k = [sourcef0; k0];

%                   nda_sourcef = .35*ones(NX-1,1);
%                    nda_sourcef = zeros(NX-1,1);
%                    nda_sourcef = sourceIteration(nda_sourcef, data);

%%%%%%%%%%%% RUN PROFUGUS SOURCE %%%%%%%%%%%%%%%%%%
%
%   Anderson = 2 norm optimization
%


if(runProfugusSource)
    
    source_mlevel = 0;
    
    prof_andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', source_mlevel);
    
    tic();
    
    [prof_sourcek, source_rhist] = anderson_knl(source_k, @profkcode, source_mlevel, prof_andoptions);
    
    toc()
    
    if(runProfugusAnderson || runProfugusNDA || runProfugusNDAAnderson || runProfugusNDA_KNL)
        s = parseXML(input_filename, running_filename, holo_elem, mat_elem);               % reset running .xml file
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%% RUN PROFUGUS ANDERSON %%%%%%%%%%%%%%%%%%

if(runProfugusAnderson)
    
    prof_andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', mlevel);
    
    tic();
    
    [prof_anders_sourcef, anders_rhist] = anderson_knl(source_k, @profkcode, mlevel, prof_andoptions);
    
    toc()
    
    if(runProfugusNDA || runProfugusNDAAnderson || runProfugusNDA_KNL)
        s = parseXML(input_filename, running_filename, holo_elem, mat_elem);               % reset running .xml file
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% RUN PROFUGUS NDA SOURCE %%%%%%%%%%%%%%%%%%

if(runProfugusNDA)
    
    source_mlevel = 0;
    
    prof_andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', source_mlevel);
    
    tic();
    
    [prof_nda_sourcek, nda_rhist] = anderson_knl(source_k, @profNDAkcode, source_mlevel, prof_andoptions);

%     [prof_nda_sourcek, nda_rhist] = anderson_knl(prof_anders_sourcef, @profNDAkcode, source_mlevel, prof_andoptions);

    toc()
    
    if(runProfugusNDAAnderson || runProfugusNDA_KNL)
        s = parseXML(input_filename, running_filename, holo_elem, mat_elem);               % reset running .xml file
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% RUN PROFUGUS NDA ANDERSON %%%%%%%%%%%%%%%%%%

if(runProfugusNDAAnderson)
    
    prof_andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', mlevel);
    
    tic();
    
    [prof_nda_anders_sourcek, nda_anders_rhist] = anderson_knl(source_k, @profNDAkcode, mlevel, prof_andoptions);
    
    toc()
    
    if(runProfugusNDA_KNL)
        s = parseXML(input_filename, running_filename, holo_elem, mat_elem);               % reset running .xml file
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


split_iter_Nps = Np * np_factor.^([0:maxit]');

iter_Nps = zeros(maxit+1, 1);

for i = 0:maxit
    if i > 0
        iter_Nps(i+1) = iter_Nps(i);
    end
    iter_Nps(i+1) = iter_Nps(i+1) + split_iter_Nps(i+1);
end

if(np_factor ~= 1)
    figure(19)
    loglog(iter_Nps(1:size(source_rhist,1)), source_rhist(:,2)/source_rhist(1,2), '-ok', ...
        iter_Nps(1:size(anders_rhist,1)), anders_rhist(:,2)/anders_rhist(1,2), '-og', ...
        iter_Nps(1:size(nda_rhist,1)), nda_rhist(:,2)/nda_rhist(1,2), '-ob', ...
        iter_Nps(1:size(nda_anders_rhist,1)), nda_anders_rhist(:,2)/nda_anders_rhist(1,2), '-or', 'LineWidth', 2)
    set(gca, 'FontSize', 14)
    xlabel(['Total Number of Histories (Np_0 = ', num2str(Np*2/1000), 'K)'])
    ylabel('Relative Residual')
    legend('Source', ['Anderson(', num2str(mlevel), ')' ], 'NDA', ['NDA-Anderson(', num2str(mlevel), ')' ])
    xlim([iter_Nps(1), iter_Nps(end)])
    ylim([1e-3, 1])
else
    figure(19)
    semilogy(iter_Nps(1:size(source_rhist,1)), source_rhist(:,2)/source_rhist(1,2), '-ok', ...
        iter_Nps(1:size(anders_rhist,1)), anders_rhist(:,2)/anders_rhist(1,2), '-og', ...
        iter_Nps(1:size(nda_rhist,1)), nda_rhist(:,2)/nda_rhist(1,2), '-ob', ...
        iter_Nps(1:size(nda_anders_rhist,1)), nda_anders_rhist(:,2)/nda_anders_rhist(1,2), '-or', 'LineWidth', 2)
    set(gca, 'FontSize', 14)
    xlabel(['Total Number of Histories (Np_0 = ', num2str(Np/1000000), 'mil)'])
    ylabel('Relative Residual')
    legend('Source', ['Anderson(', num2str(mlevel), ')' ], 'NDA', ['NDA-Anderson(', num2str(mlevel), ')' ])
    xlim([iter_Nps(1), iter_Nps(end)])
    ylim([1e-4, 1])
end

%%%%%%%%%%%% RUN PROFUGUS ANALYTIC NEWTON %%%%%%%%%%%%%%%%%%

if(runProfugusNDA_KNL)
    
    nloptions2=knl_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'x_dependent_ptv', @differential_ptv, 'p_side', 'right',...
        'x_dependent_jtv', @transport_jv, 'fx_data', 1, 'static_data', data, ...
        'etamax', -1.d-2);
    
    tic();
    
    [nda_newton_sourcef, nda_rhist, ierrk, phi_hist] = knl(sourcef0,@differential_nda_knl3,nloptions2);
    
    toc()
    
    final_res = differential_nda_knl3(nda_newton_sourcef, data);
    
    nf_res = norm(final_res, inf);
    
    nda_newton_sourcef(1:6)
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% RUN SOURCE ANDERSON %%%%%%%%%%%%%%%%%%

if(runAndersonSource)
    
    andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', mlevel);
    
    
    tic();
    
    [anders_sourcef, anders_rhist] = anderson_knl(sourcef0, @sourceIteration, mlevel, andoptions);
    
    toc()
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% RUN NDA ANDERSON %%%%%%%%%%%%%%%%%%

if(runAndersonNDA)
    
    nda_andoptions = anderson_optset('maxit',maxit,'atol',atol,'rtol', rtol,...
        'static_data', data, 'mlevel', mlevel);
    
    
    tic();
    
    [nda_anders_sourcef, nda_anders_rhist] = anderson_knl(sourcef0, @nda_sourceIteration, mlevel, nda_andoptions);
    
    toc()
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% data.ar = rec_mu;
%
% Is_zero = recoverI(nda_newton_sourcef, data, 'left')';
% Is_tau = recoverI(nda_newton_sourcef, data, 'right')';
%
% data.ar = ang_r;
%
% disp('NDA Inexact-Newton-Armijo I(0,-mu) and I(tau,mu)')
% [Is_zero(disp_mu) Is_tau(disp_mu)]
%

%
% k = size(nda_rhist,1);
% transport_sweeps = sum(nda_rhist(:,2));
%
%
%
% % if(runProfAnderson && runProfSource)
%
%
% semilogy(nda_rhist(1:k,2),nda_rhist(1:k,1)/nda_rhist(1,1),'-','Color','black','LineWidth',2);
% set(gca,'FontSize',14);
% xlabel('Transport Sweeps');
% ylabel('||r_k/r_0||');
% saveas(gcf, 'Residual_Hist','pdf')
% saveas(gcf, 'Residual_Hist','fig')
%
% plot(cgrid, nda_newton_sourcef, '-', 'Color','black','LineWidth',2)
% set(gca,'FontSize',14);
% xlabel('Space');
% ylabel('Scalar Flux');
% saveas(gcf, 'Scalar_Flux','pdf')
% saveas(gcf, 'Scalar_Flux','fig')
%
%
% plot(1:nx-1, final_res, '-', 'Color','black','LineWidth',2)
% set(gca,'FontSize',14);
% xlabel('Space');
% ylabel('Residual');
% saveas(gcf, 'Final_Residual','pdf')
% saveas(gcf, 'Final_Residual','fig')








