close all;
clear;
%% Butterworth BandPass Filter
% Courtsey desgin put on moodle by Prof. Gadre
% Removed some redundancies
%%%%%%%%
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply polynomials to get denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Band Edge speifications
fp1 = 36.8;
fs1 = 40.8;
fs2 = 60.8;
fp2 = 64.8;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lsf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lsf((B*s)/(s*s +W0*W0));     %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BSF
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BSF
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 260e3);
figure(2)
plot(f, abs(H))
title('IIR Bandstop Filter Magnitude Response')
xlabel('Frequency')
ylabel('|H(f)|')
grid 

hold on;

figure(3)
plot(f, abs(H))
title('IIR Bandstop Filter Sanity Check')
xlabel('Frequency')
ylabel('|H(f)|')
% Brickwall structure to verify design
line([0;36800],[1;1], 'Color', 'black');
line([0;36800],[0.85;0.85], 'Color', 'black');
line([40800;60800],[0.15;0.15], 'Color', 'black');
line([64800;130000],[1;1], 'Color', 'black');
line([64800;130000],[0.85;0.85], 'Color', 'black');

line([36800;36800],[0;1], 'Color', 'black');
line([64800;64800],[0;1], 'Color', 'black');
line([40800;40800],[0;0.15], 'Color', 'black');
line([60800;60800],[0;0.15], 'Color', 'black');
grid
hold off;
