function  EEG = ay_generate_sample_signal(n,rate,sec,amp_noise)
    %% Setting parameters

    % n is the number of channels
    % noise is the value of addictive noise
    % rate is the sampling rate, which means sample size per second
    % sec is the data length, which is the time
    % phs_noise is the phase noise


    %% Generate data
    EEG = zeros(n,3*sec*rate);
    dt   = 1/rate;

    %% Fully coherent with one frequency
    freq_a = 36;
    ph_a = 2*pi*rand(n,1);
    scale= 1; 
    for i = 1:n
        ind = 1:rate*sec;
        EEG(i,ind)= cos(2*pi*freq_a*dt*ind + ph_a(i))+ sqrt(scale*amp_noise)*randn(size(ind));
    end

    %% Completely non-coherent - large noise
    freq_b = 65;
    scale  = 4;
    ph_b = 2*pi*rand(n,1);
    for i  = 1:n
        ind = (1:rate*sec)+ rate*sec;
        if i<=4
            EEG(i,ind)= cos(2*pi*freq_a*dt*ind + ph_a(i))+ sqrt(scale*amp_noise)*randn(size(ind));
        else
            EEG(i,ind)= sqrt(scale*amp_noise)*randn(size(ind));
        end
    end
    
    %%  Coherent with two frequency
    ph_a  = 2*pi*rand(n,1);
    ph_b  = 2*pi*rand(n,1);
    scale = 8;
    for i = 1:n
        ind = (1:rate*sec) + 2*rate*sec;
        EEG(i,ind)= cos(2*pi*freq_a*dt*ind + ph_a(i))+cos(2*pi*freq_b*dt*ind + ph_b(i))+ sqrt(scale*amp_noise)*randn(size(ind));
    end
end 
