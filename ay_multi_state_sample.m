function Xs =ay_multi_state_sample(Ns,XSmt,SSmt,XPos,SPos,Param)
% This function samples from the psoterior estimate - smoothing result - of
% Here, we asumme that Wk is diagonal - this is part of state-transition
% model

% Model Parameters
Ak = Param.Ak;           
Wk = Param.Wk;

% Sample Result
K  = length(XSmt);
d  = size(Wk,1);
Xs = zeros(Ns,K,d);

% Find samples
for n=1:Ns
    % First sample
    %Xk = mvnrnd(XSmt{K},SSmt{K});
    Xk  = XSmt{K} + chol(SSmt{K})'*randn(d,1);
    Xs(n,K,:)=Xk'; 
    % Find a trajectory of samples    
    for k=K-1:-1:1
         Gk = SPos{k} * Ak'* pinv(Ak*SPos{k}*Ak'+Wk);
         mt = XPos{k} + Gk * ( Xk - Ak * XPos{k});
         pt = SPos{k} - Gk * (Ak*SPos{k}*Ak'+Wk)* Gk';
         %Xk = mvnrnd(mt,pt);
         Xk  = mt + chol(pt)'*randn(d,1);
         Xs(n,k,:)= Xk'; 
    end
end


