function EEG=ay_generate_signal()

ne = 4;               % number of electrodes

eig_decay_rate_a = 0.1; % higher leads to a more correlated signals
eig_decay_rate_b = 1;   % higher leads to a more correlated signals
eig_decay_rate_c = 0.1; % higher leads to a more correlated signals

coherent_freq = 59;     % one frequency
sampling_rate = 1000;   % sample per second
data_length   = 180;    % sec, will be 3 times of this

noise_var     = 0.1;
freq_shift_a  = 35;
freq_shift_b  = 0.01;
freq_shift_c  = 45;

%% generate the correlated covraince matrix
U  = orth(randn(ne,ne));

la = diag(exp(-eig_decay_rate_a*(0:ne-1)));   
Sa = U*la*U';

lb = diag(exp(-eig_decay_rate_b*(0:ne-1)));   
Sb = U*lb*U';

lc = diag(exp(-eig_decay_rate_c*(0:ne-1)));   
Sc = U*lc*U';

%% generate samples for this 300 sec
Ra = mvnrnd(zeros(ne,1),Sa);
Pa = 2*pi*rand(ne,1);

Rb = mvnrnd(zeros(ne,1),Sb);
Pb = 2*pi*rand(ne,1);

Rc = mvnrnd(zeros(ne,1),Sc);
Pc = 2*pi*rand(ne,1);

%% generate data
EEG = zeros(ne,3*data_length*sampling_rate);
for i=1:ne
    ind = 1:sampling_rate*data_length;
    dt  = 1/sampling_rate;
    EEG(i,ind)=Ra(i)*cos(2*pi*(coherent_freq+i*freq_shift_a)*dt*ind+Pa(i))+sqrt(noise_var)*randn(size(ind));
end
for i=1:ne
    ind = (1:sampling_rate*data_length)+sampling_rate*data_length;
    dt  = 1/sampling_rate;
    EEG(i,ind)=Rb(i)*cos(2*pi*(coherent_freq+i*freq_shift_b)*dt*ind+Pb(i))+sqrt(noise_var)*randn(size(ind));
end
for i=1:ne
    ind = (1:sampling_rate*data_length)+2*sampling_rate*data_length;
    dt  = 1/sampling_rate;
    EEG(i,ind)=Rc(i)*cos(2*pi*(coherent_freq+i*freq_shift_c)*dt*ind+Pc(i))+sqrt(noise_var)*randn(size(ind));
end
