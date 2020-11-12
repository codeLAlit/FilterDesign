close all;
clear;
%% Butterworth BandPass Filter
% Courtsey desgin put on moodle by Prof. Gadre
% Removed some redundancies
%%%%%%%%
Wc = 1.077;              %cut-off frequency
N = 8;                  %order 

%Poles 
p1 = Wc*cos(pi/2 + pi/16) + i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - i*Wc*sin(pi/2 + pi/16+3*pi/8);

f_samp = 330e3;         
W0 = 0.557;
B = 0.25;
[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1
%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bpf
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                         
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bpf
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 330e3);
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
% Brickwall structure to verify design
line([44000;64000],[1;1], 'Color', 'black');
line([44000;64000],[0.85;0.85], 'Color', 'black');
line([68000;165000],[0.15;0.15], 'Color', 'black');
line([0;40000],[0.15;0.15], 'Color', 'black');

line([44000;44000],[0;1], 'Color', 'black');
line([64000;64000],[0;1], 'Color', 'black');
line([68000;68000],[0;0.15], 'Color', 'black');
line([40000;40000],[0;0.15], 'Color', 'black');
grid
hold off;
