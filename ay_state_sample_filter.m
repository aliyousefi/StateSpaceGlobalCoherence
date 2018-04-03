function Xs =ay_state_sample_filter(Ns,Param,Yn)
% This function samples from the psoterior estimate - smoothing result - of
% Here, we asumme that Wk is diagonal - this is part of state-transition
% model

% Model Parameters
Ak = Param.Ak;           
Ck = Param.Ck;
Wk = Param.Wk;
X0 = Param.X0;
W0 = Param.W0;
Vk = Param.Vk;

% Sample Result
K  = length(Yn);
d  = size(Wk,1);
Xs = zeros(Ns,K,d);

% Find samples
for n=1:Ns
    % First sample
    Xk = mvnrnd(X0,W0);
    % Find a trajectory of samples    
    for k=1:K
         XPre = Ak * Xk;
         SPre = Wk;
       
         Sk   =  Ck * SPre * Ck' + Vk;
         Yp   =  Ck * XPre;
         XPos =  XPre + SPre * Ck'* Sk^-1* (Yn(k)-Yp);
         SPos = (SPre^-1 + Ck' * Vk^-1 * Ck)^-1;
               
         Xk     = mvnrnd(XPos,SPos);
         Xs(n,k,:)= Xk'; 
    end
end


