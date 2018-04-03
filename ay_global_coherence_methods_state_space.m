function [XF,SS_Param] = ay_global_coherence_methods_state_space(XF,Iter,Ns)
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


SS_Param = [];

disp('State Space')
n_channel = length(XF);
len_fs    = size(XF{1}.data,1);
% FFT Filter First
parfor p = 1:n_channel
      for j = 1:len_fs
            %disp(['Filter -ch: ' num2str(p) ' , freq: '  num2str(j)])
            %% state-space model
            Re_temp = real(XF{p}.data(j,:));%Re_temp = diff(Re_temp);
            Im_temp = imag(XF{p}.data(j,:));%Im_temp = diff(Im_temp);
            [rXSmt,rSSmt,Param,rXPos,rSPos] = ay_fft_filter([Re_temp' Im_temp'],Iter);
%            SS_Param{p}{j} = Param;
            %% draw samples
            xy     = ay_multi_state_sample(Ns,rXSmt,rSSmt,rXPos,rSPos,Param);
            Real_x = squeeze(xy(:,:,1));   
            Img_x  = squeeze(xy(:,:,2));
            for samp = 1:Ns
                XF{p}.f_data(j,samp,:) = complex(Real_x(samp,:),Img_x(samp,:));
            end
      end
end