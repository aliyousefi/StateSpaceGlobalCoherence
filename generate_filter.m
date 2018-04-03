function [filt,beta] = generate_filter(wnd_len,band_width,bins)

filt = dpss(wnd_len,band_width,1);

d_corr = xcorr(filt);
d_corr = d_corr/d_corr(length(filt));
ind    = length(filt):bins:length(d_corr);
t_corr = d_corr(ind);
beta   = exp((sum(log(t_corr)))/sum(0:length(t_corr)-1));

end