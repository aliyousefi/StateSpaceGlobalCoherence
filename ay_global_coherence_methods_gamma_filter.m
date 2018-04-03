function [COH_F,GC_Param] = ay_global_coherence_methods_gamma_filter(COH,XIter)
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

GC_Param = [];

%% GAMMA FILTER STEP
COH_F   = COH;
disp('Gamma Filter')
parfor i = 1:size(COH_F,2)
    Yn = COH(:,i);
    n  = length(Yn);
    In = ones(n,1);
    valid = ones(n,1);

    %% Set Behavioral Model and Learning Procedure
    %  Create model
    Param = ay_create_state_space(1,0,1,0,1,1,0,1,1);
    Param.Ak = 1;
    Param.Ck = 1;
    % Set learning parameters
    Param = ay_set_learning_param(Param,XIter,0,1,1,1,1,0,1,2,1);
    % EM
    [rXSmt,rSSmt,Param,rXPos,rSPos,~,EYn]=ay_em([2 0],[],In,0,Yn,[],Param,valid);

    % Extract estimating data
    COH_F(:,i)   = EYn + Param.S;
%    GC_Param{i} = Param;
end


