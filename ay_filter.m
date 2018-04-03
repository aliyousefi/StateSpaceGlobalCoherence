function [real_y,imag_y] = ay_filter(data,center_fr,lowpass_filt)
    T      = (0:length(data)-1);
    data_f = data .* exp(sqrt(-1)*2*pi*center_fr.*T);
    data_f = filter(lowpass_filt,1,data_f);
    real_y = real(data_f);
    imag_y = imag(data_f);
end