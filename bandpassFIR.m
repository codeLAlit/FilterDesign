%% FIR Bandpass filter
close all;
f_samp = 330e3; %Sampling frequency

%Band Edge speifications
fs1 = 40e3; %Stop Edge 1
fp1 = 44e3; %Pass Edge 1
fp2 = 64e3; %Pass Edge 1
fs2 = 68e3; %Stop Edge 2

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
   
%Using direct statement for Kaiser window
[n,Wn,beta,ftype] = kaiserord([fs1 fp1 fp2 fs2],[0 1 0],[0.15 0.15 0.15],f_samp);
n=n+14; %Tuning the n for desirable results
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); % impulse response in time domain
kaiser_coeffs=kaiser(n+1, beta); %Kaiser window coefficients
fvtool(hh); %Analysing the Frequency and Phase response



%Magnitude Response
[H,f] = freqz(hh,1,1024, f_samp); % Magnitude respnse
figure(2)
title("FIR Bandpass Filter Magnitude Response")
xlabel("Freqency")
ylabel("|H(f)|")
plot(f,abs(H))
grid
hold on;
figure(3)
plot(f, abs(H))
title("FIR Bandpass Filter Sanity Check")
xlabel("Freqency")
ylabel("|H(f)|")
%Creating Brick wall diagram to verify design
line([44000;64000],[1.15;1.15], 'Color', 'black');
line([44000;64000],[1;1], 'Color', 'black');
line([44000;64000],[0.85;0.85], 'Color', 'black');
line([68000;165000],[0.15;0.15], 'Color', 'black');
line([0;40000],[0.15;0.15], 'Color', 'black');

line([44000;44000],[0;1.15], 'Color', 'black');
line([64000;64000],[0;1.15], 'Color', 'black');
line([68000;68000],[0;0.15], 'Color', 'black');
line([40000;40000],[0;0.15], 'Color', 'black');
grid

hold off;