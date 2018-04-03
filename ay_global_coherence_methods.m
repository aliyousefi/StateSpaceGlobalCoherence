function [COH,COH_F,SS_Param,GC_Param] = ay_global_coherence_methods(method,EEG,Wnd,Iter,fs,Ns,XIter,bins)
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
SS_Param = [];

%% generate slepian window
base_filt=generate_filter(Wnd,2.5,Wnd);

%% run FFT
disp('Filter')
n_channel = size(EEG,1);
XF        = cell(n_channel,1);
parfor c  = 1:n_channel
    for j = 1:length(fs)
        [temp_a,temp_b]   = ay_filter(EEG(c,:),fs(j),base_filt);
        XF{c}.data(j,:)   = temp_a(1:Wnd:end)+sqrt(-1)*temp_b(1:Wnd:end);
    end
end


Ms = floor(length(XF{1}.data(1,:)));
%% PNAS method - there is no de-mean
if method == 1 || method == 0
    COH = zeros(floor(Ms/bins),length(fs));
    ms  = bins;
    ind_fill = 0; 
    for i = 1:ms:Ms-ms
      ind_fill =  ind_fill +1;
      ind_a = i;
      ind_b = i+ms-1 ;
      for j = 1:length(fs)
            cov = zeros(n_channel,n_channel);
            for p=1:n_channel
               Ta = XF{p}.data(j,ind_a:ind_b)';
               if method==1
                rTa = real(Ta)-mean(real(Ta));
                iTa = imag(Ta)-mean(imag(Ta));
                Ta  = rTa + iTa*sqrt(-1);
               end
             
                for q=p:n_channel
                    Tb = XF{q}.data(j,ind_a:ind_b)';
                    if method==1
                      rTb = real(Tb)-mean(real(Tb));
                      iTb = imag(Tb)-mean(imag(Tb));
                      Tb  = rTb + iTb*sqrt(-1);
                    end
                    cov(p,q) = sum(Ta.*conj(Tb));
                    cov(q,p) = conj(cov(p,q));
                end
            end
        eigs    = svd(cov);
        COH(ind_fill,j)= eigs(1)/sum(eigs);
      end
    end
else
    %% COH matrix
    COH  = zeros(size(XF{1},1),length(fs));
    ms   = 1;
    %% State Space Method
    %  run this for each channel and different frequencies 
    disp('State Space')
    % FFT Filter First
    parfor p = 1:n_channel
      for j = 1:length(fs)
            %disp(['Filter -ch: ' num2str(p) ' , freq: '  num2str(j)])
            %% state-space model
            Re_temp = real(XF{p}.data(j,:));
            Im_temp = imag(XF{p}.data(j,:));
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
    
    %% Global Coherence Measure
    % we have multiple measure of FFT per Wnd - number: Wnd/Ls
    disp('Covariance Step')
    
    %% only scaled
    if method ==2
        for i = 1:Ms
         %    disp(['COV:' num2str(i) ' out of ' num2str(Ms) ])  
            ind_a  = (i-1)*ms+1;
            ind_b  = i*ms ;
            
            for j = 1:length(fs)
                cov = zeros(n_channel,n_channel);
                for p=1:n_channel
                    Ta   = squeeze(XF{p}.f_data(j,:,ind_a:ind_b))';
                    % method 1:take one sample and center by the mean of phase
                    da   = sqrt(sum(Ta.*conj(Ta)));
                    for q=p:n_channel
                        Tb   = squeeze(XF{q}.f_data(j,:,ind_a:ind_b))';
                        % method 1:take one sample and center by the mean of phase
                        db   = sqrt(sum(Tb.*conj(Tb)));
                        cov(p,q) = sum(Ta.*conj(Tb))/(da*db);
                        cov(q,p) = conj(cov(p,q));
                    end
                end
                eigs    = svd(cov);
                COH(i,j)= eigs(1)/sum(eigs);
            end
        end
    end
    
    %% nothing
    if method ==3
        for i = 1:Ms
         %    disp(['COV:' num2str(i) ' out of ' num2str(Ms) ])  
            ind_a  = (i-1)*ms+1;
            ind_b  = i*ms ;
            
            for j = 1:length(fs)
                cov = zeros(n_channel,n_channel);
                for p=1:n_channel
                    Ta   = squeeze(XF{p}.f_data(j,:,ind_a:ind_b))';
                    % method 1:take one sample and center by the mean of phase
                    for q=p:n_channel
                        Tb   = squeeze(XF{q}.f_data(j,:,ind_a:ind_b))';
                        % method 1:take one sample and center by the mean of phase
                        cov(p,q) = sum(Ta.*conj(Tb));
                        cov(q,p) = conj(cov(p,q));
                    end
                end
                eigs    = svd(cov);
                COH(i,j)= eigs(1)/sum(eigs);
            end
        end
    end
    
    %% de-meaned, scaled
    if method == 4 
        for i = 1:Ms
         %    disp(['COV:' num2str(i) ' out of ' num2str(Ms) ])  
            ind_a  = (i-1)*ms+1;
            ind_b  = i*ms ;
            
            for j = 1:length(fs)
                cov = zeros(n_channel,n_channel);
                for p=1:n_channel
                    Ta   = squeeze(XF{p}.f_data(j,:,ind_a:ind_b))';
                    % method 1:take one sample and center by the mean of phase
                    Ta = (real(Ta)-mean(real(Ta)))+(imag(Ta)-mean(imag(Ta)))*sqrt(-1);
                    da   = sqrt(sum(Ta.*conj(Ta)));
                    for q=p:n_channel
                        Tb   = squeeze(XF{q}.f_data(j,:,ind_a:ind_b))';
                        % method 1:take one sample and center by the mean of phase
                        Tb = (real(Tb)-mean(real(Tb)))+(imag(Tb)-mean(imag(Tb)))*sqrt(-1);
                        db   = sqrt(sum(Tb.*conj(Tb)));
                   
                        cov(p,q) = sum(Ta.*conj(Tb))/(da*db);
                        cov(q,p) = conj(cov(p,q));
                    end
                end
                eigs    = svd(cov);
                COH(i,j)= eigs(1)/sum(eigs);
            end
        end
    end
    
    %% de-meaned, not scaled
    if method == 5
        for i = 1:Ms
         %    disp(['COV:' num2str(i) ' out of ' num2str(Ms) ])  
            ind_a  = (i-1)*ms+1;
            ind_b  = i*ms ;
            
            for j = 1:length(fs)
                cov = zeros(n_channel,n_channel);
                for p=1:n_channel
                    Ta   = squeeze(XF{p}.f_data(j,:,ind_a:ind_b))';
                    % method 1:take one sample and center by the mean of phase
                    Ta = (real(Ta)-mean(real(Ta)))+(imag(Ta)-mean(imag(Ta)))*sqrt(-1);
                   
                    for q=p:n_channel
                        Tb   = squeeze(XF{q}.f_data(j,:,ind_a:ind_b))';
                        % method 1:take one sample and center by the mean of phase
                        Tb = (real(Tb)-mean(real(Tb)))+(imag(Tb)-mean(imag(Tb)))*sqrt(-1);
                        cov(p,q) = sum(Ta.*conj(Tb));
                        cov(q,p) = conj(cov(p,q));
                    end
                end
                eigs    = svd(cov);
                COH(i,j)= eigs(1)/sum(eigs);
            end
        end
    end

end

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


