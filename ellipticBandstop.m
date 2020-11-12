%% Elliptical Bandstop filter
clear 
close all;
%Low pass transformed parameters
fp = 1; fs = 1.367; D1=0.384; D2=43.44;

%Calculating important parameters
Wp = 2*pi*fp; 
Ws = 2*pi*fs;
ep = sqrt(D1); 
es = sqrt(D2);
rp=10*log10(1+D1);
rs=10*log10(1+D2);
k = Wp/Ws;
k1 = ep/es;

% Elliptic Integral calculations
[K,Kp] = ellipk(k);
[K1,K1p] = ellipk(k1);

%Order Calcaulation
Nexact = (K1p/K1)/(Kp/K);  
N = ceil(Nexact);

%Poles and Zeros of H(sl)
[zeros, poles, dcgain]=ellipap(N, rp, rs);
[num,den] = zp2tf(zeros, poles, dcgain);

f_samp = 260e3;  %sampling frequency       
W0 = 0.688; % W0 frequency as defined
B = 0.52;   % B parameter as defined

%Defining symbolic polynomial so that all the transformations are
% followed like the theory
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bpf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bpf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure(2)
plot(f, abs(H))
title('IIR Bandpass Filter Magnitude Response')
xlabel('Frequency')
ylabel('|H(f)|')
grid 

hold on;

figure(3)
plot(f, abs(H))
title('IIR Bandpass Filter Sanity Check')
xlabel('Frequency')
ylabel('|H(f)|')
%Brick wall structure to verify design
line([0;36800],[1;1], 'Color', 'black');
line([0;36800],[0.85;0.85], 'Color', 'black');
line([40800;60800],[0.15;0.15], 'Color', 'black');
line([64800;130000],[1;1], 'Color', 'black');
line([64800;130000],[0.85;0.85], 'Color', 'black');

line([36800;36800],[0;1], 'Color', 'black');
line([64800;64800],[0;1], 'Color', 'black');
line([40800;40800],[0;0.15], 'Color', 'black');
line([60800;60800],[0;0.15], 'Color', 'black')
grid
hold off;
