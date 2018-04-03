function [COH] = ay_global_coherence_methods_global_coherence(method,XF,bins)
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

n_channel = length(XF);
%% PNAS method - there is no de-mean
if method == 1 || method == 0
    len_fs  = size(XF{1}.data,1);
    Ms      = floor(length(XF{1}.data(1,:)));

    COH = zeros(floor(Ms/bins),len_fs);
    ms  = bins;
    ind_fill = 0; 
    for i = 1:ms:Ms-ms
      ind_fill =  ind_fill +1;
      ind_a = i;
      ind_b = i+ms-1 ;
      for j = 1:len_fs
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
    len_fs  = length(squeeze(XF{1}.f_data(:,1,1)));
    Ms      = floor(length(squeeze(XF{1}.f_data(1,1,:))));

    %% COH matrix
    COH  = zeros(Ms,len_fs);
    ms   = 1;
   
    
    %% Global Coherence Measure
    % we have multiple measure of FFT per Wnd - number: Wnd/Ls
    disp('Covariance Step')
    
    %% only scaled
    if method ==2
        for i = 1:Ms
         %    disp(['COV:' num2str(i) ' out of ' num2str(Ms) ])  
            ind_a  = (i-1)*ms+1;
            ind_b  = i*ms ;
            
            for j = 1:len_fs
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
            
            for j = 1:len_fs
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
            
            for j = 1:len_fs
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
            
            for j = 1:len_fs
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

