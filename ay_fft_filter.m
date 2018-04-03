function [XSmt,SSmt,Param,XPos,SPos]= ay_fft_filter(Yn,Iter)

%% Parameters
Ak  = eye(2,2);
Ck  = eye(2,2);
%% Parameters
Param.Wk = cov(diff(Yn)).*eye(2,2);
Param.Wk = 1e-3*0.5*(Param.Wk+Param.Wk'); 
Param.Vk = cov(diff(Yn)).*eye(2,2);
Param.Vk = 0.5*(Param.Vk+Param.Vk');
Param.W0 = 1e2*Param.Wk;
Param.X0 = [0;
            0];

 
K = size(Yn,1);
for iter=1:Iter
    %% data length
    K = size(Yn,1);
    %% variables
    % one-step mean and variance
    XPre = cell(K,1);
    SPre = cell(K,1);
    % filter mean and covariance
    XPos = cell(K,1);
    SPos = cell(K,1);
    % As
    As   = cell(K,1);
    % smoother mean and covariance
    XSmt = cell(K,1);
    SSmt = cell(K,1);
    % extra components
    Ckk   = cell(K,1);
    Ckk_1    = cell(K,1);
    Wkk_1    = cell(K,1);
    
    %% Parameters
    Wk = Param.Wk;
    Vk = Param.Vk;
    W0 = Param.W0;
    X0 = Param.X0;

    
    %% run filter step
    for k=1:K
        %% Run one step prediction
        if k == 1
            XPre{k} = Ak * X0 ;
            SPre{k} = Ak * W0 * Ak'+ Wk;
        else
             XPre{k} = Ak * XPos{k-1};
             SPre{k} = Ak * SPos{k-1}* Ak' + Wk;
        end
        %% Run filter
        Yp =  Yn(k,:)' - Ck * XPre{k};
        Sk =  Ck * SPre{k} * Ck' + Vk;
        Kk =  SPre{k}* Ck' * Sk^-1;
        XPos{k} =  XPre{k} + Kk * Yp;
        SPos{k} = (eye(2,2)-Kk*Ck)*SPre{k};
%         
%          Sk      =  Ck * SPre{k} * Ck' + Vk;
%          Yp      =  Ck * XPre{k} ;
%          XPos{k} =  XPre{k} + SPre{k} * Ck'* pinv(Sk)* (Yn(k,:)'-Yp);
%                     % SPos
%          SPos{k} = pinv(pinv(SPre{k}) + Ck' * pinv(Vk) * Ck);
    end
    
    %% run smoother step
    XSmt{end} = XPos{end};
    SSmt{end} = SPos{end};
    for k=K-1:-1:1
        % Ak, equation (A.10)
        As{k} = SPos{k} * Ak' *  SPre{k+1}^-1 ;
        % Smting function, equation (A.9)
        XSmt{k} = XPos{k} + As{k} * (XSmt{k+1}- XPre{k+1});
        % Variance update, equation (A.11)
        SSmt{k} = SPos{k} + As{k} * (SSmt{k+1}- SPre{k+1}) * As{k}';
    end
    % Kalman smoother for time 0
    As0   = W0 * Ak' *  SPre{1}^-1;
    XSmt0 = X0 + As0 * (XSmt{1}- XPre{1}) ;
    SSmt0 = W0 + As0 * (SSmt{1}- SPre{1}) * As0';
    
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk)
    % Ckk = E(Xk*Xk) prediction by smoothing 
    for k = 1:K
        % Wk update - Smting Xk*Xk
        Ckk{k} = SSmt{k} + XSmt{k} * XSmt{k}';
    end
    Ckk0 = SSmt0 + XSmt0 * XSmt0';
    
    %% Extra Component of the State Prediction Ckk = E(Xk*Xk-1)
    % Ckk_1=E(Xk-1*Xk) prediction by smoothing - it is kept at index K
    % Wkk_1= Ckk_1 + Bias
    % Covariance for smoothed estimates in state space models - Biometrica
    % 1988- 601-602
    for k = 1:K
        % Wkk update - Smoothing Xk-1*Xk
        if k>1
            Wkk_1{k}   = As{k-1} * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt{k-1}'; 
        else
            Wkk_1{k}   = As0 * SSmt{k};
            Ckk_1{k}   = Wkk_1{k} + XSmt{k} * XSmt0'; 
        end
    end
    
    %% Parameter Update
    Wk = Ckk{1}+Ckk0-(Ckk_1{1}+Ckk_1{1}');
    % state process parameters
     for k = 2:K
        Wk =  Wk+ Ckk{k}+Ckk{k-1}-(Ckk_1{k}+Ckk_1{k}');
     end
     
    Wk = (Wk+Wk')/(2*K);
   % Wk = Wk/K;
    Wk = Wk .*eye(2,2);
    
    % observation process is independent
    Vk = zeros(2,2);
    for k = 1:K
       Vk = Vk + (Yn(k,:)'-XSmt{k})*(Yn(k,:)'-XSmt{k})' + SSmt{k};
    end
    Vk = (Vk+Vk')/(2*K);
    Vk = Vk.*eye(2,2);

    % set X0 and W0
    X0 = XSmt0;
    W0 = K*Wk;
    
    %% Parameters
    Param.Wk = Wk;
    Param.Vk = Vk;
    Param.W0 = W0;
    Param.X0 = X0;
end        
Param.Ak = Ak;
Param.Ck = Ck;

