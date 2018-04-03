function XF = ay_fft(method,X,Ls,Wnd)
% Calculate FFT for each Wnd in chunk of Ls
% Here, there is no smoothing function
% method 1 is un-normalized energ
% method 2 is normalized energy

if method ==1
    U = window('hann',Ls); U = U/sqrt(sum(U.^2));
else
    U = ones(Ls,1);
end

% XF row is number of Wnd segments, col is number of Ls segments in Wnd
XF    = cell(floor(length(X)/Wnd),floor(Wnd/Ls));
[m,n] = size(XF);
for i=1:m
    for j=1:n
        ind_a  = (i-1)*Wnd+(j-1)*Ls+1;
        ind_b  = ind_a+Ls-1;
        tX     = ay_detrend(X(ind_a:ind_b)');
        tX     = tX.*U;
%         if method > 2    % only in PLV method
%                 tX=tX/sqrt(sum(tX.^2));
%         end
         XF{i,j}= fft(tX)/length(tX); 
    end
end