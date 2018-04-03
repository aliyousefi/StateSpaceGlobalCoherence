function [XF,base_filt] = ay_global_coherence_methods_filter(EEG,Wnd,fs)
%% Input Argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method : 1 is PNAS paper with Gamma Filter (no energy normalization, hamming window)
%        : 2 based on state-space and sampling method (scaled)
%        : 3 based on state-space and sampling method (nothing)
%        : 4 based on state-space and sampling method (de-meaned, scaled)
%        : 5 based on state-space and sampling method (de-meaned)

% bins   : should 1 for methods 2,3,4,5; for method 1, it defines number of windows

% Ns     : number of samples used for used for cross spectral matrix est. (method 2,3,4,5)

% fs     : global coehrency is measured at these frequencies (int/Wnd)

% Iter   : iteration of EM for state space model

% XIter  : iteration for gamma filter EM

% Wnd    : length of FFT window

% EEG    : signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output Argument
% COH    : coherence measure
% COH_F  : filtered coherence
% SS_Param   : state space model param    
% GC_Param  : global coherence filter param


%% generate slepian window
base_filt=generate_filter(Wnd,2.5,Wnd);

%% run FFT
disp('Filter')
n_channel = size(EEG,1);
XF        = cell(n_channel,1);
parfor c  = 1:n_channel
    for j = 1:length(fs)
        [temp_a,temp_b]   = ay_filter(EEG(c,:),fs(j),base_filt);
        XF{c}.data(j,:)   = temp_a(Wnd:Wnd:end)+sqrt(-1)*temp_b(Wnd:Wnd:end);
    end
end
