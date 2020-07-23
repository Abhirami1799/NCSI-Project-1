%%  Question-1

%% Initializing some terms

BF = [500, 4000];
intensities = (-10:10:80);
h = (0:1/8:9);
tones = zeros(1,73);
for i = 1:73
    tones(i) = 62.5*(2.0^h(i));
end

rt = 10e-3; % ramp time in seconds
T = 200e-3; % stimulus duration in seconds
Fs = 100e3; % sampling rate in Hz (must be 100, 200 or 500 kHz)
% PSTH parameters
nrep = 50;  % number of stimulus repetitions
psthbinwidth = 0.5e-3; % binwidth in seconds
irpts = rt*Fs;
t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);

rate_matrix = zeros(73,10,2); % (frequency, intensity, BF)


%% Creating the stimuli and generating output

for i = 1:2
% model fiber parameters
CF    = BF(uint8(i))   % CF in Hz;   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
    for j = 1:10
    % stimulus parameters
    stimdb = intensities(j);  % stimulus intensity in dB SPL
        for k = 1:73
            
            disp("Int - " + j + ". Freq - " + k);
            
            F0 = tones(k);     % stimulus frequency in Hz   

            pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
            pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
            pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

            vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
            [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
            
            rate_matrix(k,j,i) = sum( psth( 1: (length(psth)/2) ) ); 
            timePSTH = length(psth)/Fs;
            
        end
    end
end
rate_matrix = rate_matrix / nrep;
rate_matrix = rate_matrix / ((length(psth)/Fs)/2);


%% Plotting Rate vs Frequency graphs

intensities = (-10:10:80);
tonelog = log10(tones);
for i = 1:2
    figure(i);
    grid on;
    hold on;
    
    for j = 1:10
        plot(tonelog, squeeze( rate_matrix(:,j,i) ) );
    end
    xlabel('Frequency');
    ylabel('Rate');
    title(['Rate vs Frequency ', num2str(BF(i)), 'Hz ANF']);
    hold off;
end


%% Plotting Rate vs Intensity curves

rate_int = zeros(10,2);
for i = 1:2
    for j = 1:10
        rate_int(j,i) = rate_matrix(25*i-(i-1),j,i);
    end
end
figure(3)
hold on;
grid on;
plot((-10:10:80),rate_int(:,1));
plot((-10:10:80),rate_int(:,2));
xlabel('Intensity');
ylabel('Rate');
title('Rate vs Intensity for both ANFs');
legend('500 Hz','4000 Hz');


%% Question - 2

%% Separating out 'ah', generating different intensity stimuli 

h = (0:1/8:7);
ANF = zeros(1,length(h));
for k = 1:length(h)
    ANF(k) = 62.5*(2^h(k));
end

intensities = -20:10:80;
[fivewo, fs] = audioread('fivewo.wav');
fivewo = fivewo' ;
S = 1+round(1.05*fs);
E = 1+round(1.15*fs);
L = length(fivewo);
[y,fs] = audioread('fivewo.wav',[1+round(1.05*fs),1+round(1.15*fs)]);
y = y';
Fs = 100e3;
rt = 10e-3;
%sound(y,Fs);
tester = [fivewo(1:S) zeros(1,50000) fivewo(S:E) zeros(1,50000) fivewo(E:L)] ;
sound(tester,Fs);
%t = 0:1/Fs:0.2*fs/Fs; % time vector
t = (0:length(y)-1)/Fs;
mxpts = length(t);
irpts = rt*Fs;
y(1:irpts) = y(1:irpts) .* ((0:(irpts-1))/irpts); 
y((mxpts-irpts):mxpts) = y((mxpts-irpts):mxpts).*((irpts:-1:0)/irpts);
%s = size(t);
x = rms(y);
I_rms = 20*log10(x/(20*10^(-6)));
rate_matrix2 = zeros(1,11);
Input = zeros(1,length(y));
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3;
implnt = 0;
nrep = 50;
CF = 500;

psthbinwidth = 0.5e-3;
%psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
ahTime = length(y) / Fs ;
for i = 1:11 
    Input = y * 10^(((i-1)*10 - 95)/20);
    vihc = catmodel_IHC(Input,CF,nrep,1/Fs,ahTime*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
    %{
    timeout = (1:length(psth))*1/Fs;
    psthtime = timeout(1:psthbins:end); % time vector for psth
    pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
    Psth = pr/psthbinwidth; % psth in units of spikes/s
    %}
    rate_matrix2(i) = sum(psth(1:(length(psth)/2)));
end
rate_matrix2 = rate_matrix2 / nrep;
rate_matrix2 = rate_matrix2 / ahTime ;


%% Plotting Rate vs Intensity of 'ah' and that of 500 Hz ANF for stimulus at BF

figure(4)
grid on;
hold on;
plot((-20:10:80), rate_matrix2);         % for the ah stimuylus
plot((-10:10:80), rate_matrix(25,:,1));  % for the 500Hz ANF, stim @ BF
xlabel('Intensity');
ylabel('Rate');
title('Rate vs Intensity');
legend('ah','500Hz');


%% Recording for all ANFs at 3 dB levels

% taking the 3 intensities as 0dB, 35dB, and 70dB
fwTime = length(fivewo)/Fs;
rate_0 = zeros(54,fwTime*Fs*2);
rate_35 = zeros(54,fwTime*Fs*2);
rate_70 = zeros(54,fwTime*Fs*2);

% playing fivewo for the ANFs
for i = 4:57
    CF = ANF(i);
    Input = fivewo .* 10^(-75/20);
    vihc = catmodel_IHC(Input,CF,nrep,1/Fs,fwTime*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
    rate_0(i-3,:) = psth;
end

for i = 4:57
    CF = ANF(i);
    Input = fivewo .* 10^(-45/20);
    vihc = catmodel_IHC(Input,CF,nrep,1/Fs,fwTime*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
    rate_35(i-3,:) = psth;
end

for i = 4:57
    CF = ANF(i);
    Input = fivewo .* 10^(-5/20);
    vihc = catmodel_IHC(Input,CF,nrep,1/Fs,fwTime*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
    rate_70(i-3,:) = psth;
end


%% Plotting the actual spectrogram and the average ANF response rates

figure(5)
spectrogram(fivewo, hann(25.6e-3*Fs), 12.8e-3*Fs, 1:8000, Fs, 'yaxis');
title('Spectrogram for the speech signal');

p = (2:7);
wind = zeros(1,6);
for i = 1:length(p)
    wind(i) = 1e-3 * Fs * (2^p(i));
end

% wind(i) has the number of samples in each window
winShift = floor(wind/2);

F = ANF(4:57);
for w = 1:6      % for each window size
    t2 = wind(w)/2 : winShift(w) : length(fivewo)-wind(w)/2;
    avg_rates = zeros(length(F), length(t2));
    for f = 1:54 % for each ANF CF
        for i = 1:length(t2) % b for bin number
            xo = rate_70(f,(t2(i)-winShift(w)+1):(t2(i)+winShift(w)));
            avg_rates(f,i) = sum(xo)*Fs / wind(w);
        end
    end
    figure(6);
    subplot(2,3,w);
    
    [ tim, frq ] = meshgrid( t2, F);
    surf(tim, frq, avg_rates/nrep,'edgecolor','none');
    set(gca,'xtick',[]);set(gca,'ytick',[]);xlabel([]);ylabel([]);colorbar('off');
    xlim([0,1.5e5]);
    title(['Window Size = ',num2str(wind(w)/(1e-3*Fs)),'ms']);
    xlabel('Time');
    ylabel('Frequency');
    view(2);
end


%% Question - 3

%% Generating FFTs for small time windows, finding the frequency corresponding to max energy, plotting on spectrogram

figure(7)
spectrogram(fivewo, hann(12.8e-3*Fs), 6.4e-3*Fs, 1:8000, Fs, 'yaxis');
view(3);
hold on;

nF =  [6, 10, 14, 18, 22, 26, 30, 34, 38];
cmap1 = hsv(9);
win = 12.8e-3*Fs;
wshift = floor(win/2);
t3 = win/2 : wshift : length(fivewo)-win/2;
fre_pt = zeros(1,length(t3));

for f = 1 : 9
    for i = 6 : length(t3)
        Xp = rate_70(nF(f), (t3(i)-wshift+1) : (t3(i)+wshift));
        m = mean(Xp);
        FFT = abs(fft(Xp - m));
        [M,I] = max(squeeze(FFT(1:length(FFT)/2)));
        fre_pt(i) = I*Fs/length(FFT);
    end
    %{
    line = zeros(1,length(t3));
    for k = 1 : length(t3)
        line(k) = ANF(nF(f));
    end
    plot(t3,line);
    hold on;
    %}
    figure(7);
    scatter3(t3/Fs,fre_pt/1000,zeros(1,length(t3))-10,[],cmap1(f,:),'filled', 'MarkerEdgeColor', 'k');
    ylim([0 3]);
    hold on;
end
view(2);