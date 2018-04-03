clear all
% %% Generate data BY Dr code
load('eeganes07laplac250_detrend_all.mat')
sample_rate = 250;     % Sampling Rate 
data_length = 9267.8;  % Data Length in Seconds


%% change to 1 for all channels
ch_jump    = 32; % 1;

%% Model Setting for PNAS
% change to a larger number - 25 or 100 - to guarantee EM will converge
StateIter  = 5;   %25;
% change to a larger number - 100 - to guarantee filter will converge
GammaIter  = 10;  %100;
% change to a larger number - twice as large as number of channels - to have an unbiased covariance matrix
Ns         = 32; %128;
% FFT Wnd, you can set the window to 1024
Wnd        = 2048;
% Set frequency accordingly
Freq       = (10:150)/Wnd;  
% You can define which part of the data to be processed
part_data  = [2.5e4   1e5]; %[2.5e4   size(data_det,1)-2.5e4];
kk =1;
data_det_x = data_det(part_data(kk,1):part_data(kk,2),1:ch_jump:end);

%% Run Filter First
[XF,base_filt]  = ay_global_coherence_methods_filter(data_det_x',Wnd,Freq);

%% Method 0 - PNAS
method = 0;
Bins   = 8; 
COH = ay_global_coherence_methods_global_coherence(method,XF,Bins);
[COH_F,ParamFilter] = ay_global_coherence_methods_gamma_filter(COH,GammaIter);
save(['Part_0_' num2str(kk)],'data_det','data_det_x','COH','COH_F');

%% PNAS without mean dropped
method = 1;
COH    = ay_global_coherence_methods_global_coherence(method,XF,Bins);
[COH_F,ParamFilter] = ay_global_coherence_methods_gamma_filter(COH,GammaIter);
save(['Part_1_' num2str(kk)],'COH','COH_F');


%% State Space Call
[XF,ParamState] = ay_global_coherence_methods_state_space(XF,StateIter,Ns);


%% Model Setting for State Space Model Method 2
method = 2;
Bins   = 1;
COH    = ay_global_coherence_methods_global_coherence(method,XF,Bins);
[COH_F,ParamFilter] = ay_global_coherence_methods_gamma_filter(COH,GammaIter);
save(['Part_2_' num2str(kk)],'COH','COH_F');

%% Model Setting for State Space Model Method 3
method = 3;
COH    = ay_global_coherence_methods_global_coherence(method,XF,Bins);
[COH_F,ParamFilter] = ay_global_coherence_methods_gamma_filter(COH,GammaIter);
save(['Part_3_' num2str(kk)],'COH','COH_F');
