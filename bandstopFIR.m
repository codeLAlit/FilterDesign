%% FIR Bandpass filter
close all;
f_samp = 260e3; %Sampling frequency

%Band Edge speifications
fs1 = 40.8e3; %Stop Edge 1
fp1 = 36.8e3; %Pass Edge 1
fp2 = 64.8e3; %Pass Edge 2
fs2 = 60.8e3; %Stop Edge 2

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
   
%Using direct statement for Kaiser window
[n,Wn,beta,ftype] = kaiserord([fp1 fs1 fs2 fp2],[1 0 1],[0.15 0.15 0.15],f_samp);
n=n+14; %Tuning the n for desirable results
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); % impulse response in time domain 
kaiser_coeffs=kaiser(n+1, beta); %Kaiser window coefficients
fvtool(hh); %Analysing the Frequency and Phase response



%magnitude response
[H,f] = freqz(hh,1,1024, f_samp); % Magnitude respnse
figure(2)
title("FIR Bandstop Filter Magnitude Response")
xlabel("Frequency")
ylabel("|H(f)|")
plot(f,abs(H))
grid
hold on;

figure(3)
plot(f, abs(H))
title("FIR Bandstop Filter Sanity Check")
xlabel("Frequency")
ylabel("|H(f)|")
%Creating Brick wall diagram to verify design 
line([0;36800],[1.15;1.15], 'Color', 'black');
line([0;36800],[1;1], 'Color', 'black');
line([0;36800],[0.85;0.85], 'Color', 'black');
line([40800;60800],[0.15;0.15], 'Color', 'black');
line([64800;130000],[1;1], 'Color', 'black');
line([64800;130000],[0.85;0.85], 'Color', 'black');
line([64800;130000],[1.15;1.15], 'Color', 'black');

line([36800;36800],[0;1.15], 'Color', 'black');
line([64800;64800],[0;1.15], 'Color', 'black');
line([40800;40800],[0;0.15], 'Color', 'black');
line([60800;60800],[0;0.15], 'Color', 'black');
grid

hold off;