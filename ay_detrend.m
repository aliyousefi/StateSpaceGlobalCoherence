function Xd =ay_detrend(X)
% this function detrends signal
p  = polyfit([1:length(X)]',X,1);
Xd = X-polyval(p,[1:length(X)]');