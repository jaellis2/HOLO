function [phiHO, ephiHO, jHO] = running_HO(phi, data, Np_update, newfilename)


%FIX FOLDER OPERATIONS, CURRENTLY HARD CODED FOR NCSU'S HPC CLUSTER


hdf5file = [data.prob_name,'_flux.h5'];
j_hdf5file = [data.prob_name,'_current.h5'];




h = data.tau/(length(phi));

lin_source = ((data.csig_s.*phi) + data.q); % LINEAR NEUTRON DENSITY

% vol_source = lin_source*h^2;

vol_source = lin_source;

phi_norm = (h^3)*sum(vol_source); % HOW TO CALCULATE
% normal_vol_source = vol_source / phi_norm; % PROFUGUS SCALES AS WELL

worked = addSource(newfilename, vol_source, Np_update, [0;0], data);


%%% TO CHANGE FOR PC's, RIGHT NOW ONLY
%%% FOR NCSU HPC


%%%%%%%%%%%%%%%%
%  RUN BSUB
%




% s = system('bsub < /home/jaellis2/CASL/HOLO/runProfMC')    %NCSU HPC

runcommand = [data.xmc_pwd,'xmc -i ', data.newfilename];

% s = system('~/Documents/SOFTWARE/ProfugusB/bin/xmc -i running_mcfixedsrc_cart.xml') % Testing purposes

s = system(runcommand)

if(s ~= 0)
    return;
end

% outputfile = '.hdf5';
% 
% k = 1;
% while(~exist(hdf5file, 'file'))
% 
% pause(5)
% 
% k = k+1;
% 
% if(k > 600)
%    quit;
% end
% end

[uns_vol_phiHO, uns_face_ephiHO, uns_face_jHO] = extractSource(hdf5file, j_hdf5file, data);


% s = system('mv ~/Documents/SOFTWARE/HOLO/*.h5 ~/Documents/SOFTWARE/HOLO/HighOrderSolutions/');

mvcommand = ['mv *.h5 ', data.HO_sol_pwd];
s = system(mvcommand);

if(s ~= 0)
    return; 
end

%  phiHO = h*uns_phiHO;
% ephiHO = h*uns_ephiHO;
% jHO = h*uns_jHO;                     % QUESTIONNNNNNN
% 
% phiHO   = phi_norm / h^2 * uns_vol_phiHO;
% ephiHO = phi_norm / h * uns_vol_ephiHO;
% jHO       = phi_norm / h * uns_vol_jHO;                     % QUESTIONNNNNNN



vol_phiHO = phi_norm * uns_vol_phiHO;
face_ephiHO = phi_norm * uns_face_ephiHO;
face_jHO = phi_norm * uns_face_jHO;

% phiHO = vol_phiHO*h^2;
% ephiHO = face_ephiHO*h;
% jHO = face_jHO*h;

phiHO = vol_phiHO;
ephiHO = face_ephiHO;
jHO = face_jHO;


%

% worked = addSource(newfilename, filename, sourcef);

figure(2)
plot([.01:.02:.99], phiHO, 'g-', [0:.02:1], ephiHO, 'b-')
legend('Cell', 'Edge')

figure(3)
plot( [0:.02:1], jHO, 'm-')
legend('Current')

% avg = (phiHO + (ephiHO(1:end-1) + ephiHO(2:end))/2) /2;
% 
% plot(avg)
% legend('Cell, Edge Averaged')
% 
% Test = open('Test.fig');
%     figure(Test)
%      pause(2)  %pause for you to see that it is bringing ur fig up
%     hold on
%     plot(avg)
% saveas(gcf,'Test', 'fig')

end