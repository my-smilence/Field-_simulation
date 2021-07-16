clc
clear
%% %%%%%%%%%%%%%%%%%%%%     
%  Set initial parameters
f0=3e6; %  Transducer center frequency [Hz]
fs=100e6; %  Sampling frequency [Hz]
c=1540; %  Speed of sound [m/s]
lambda=c/f0; %  Wavelength [m]
height=5/1000; %  Height of element [m]
width=1/1000; %  Width of element [m]
kerf=width/5; %  Distance between transducer elements [m]
N_elements=3; %  Number of elements
N_elements2=16; %  Number of elements
focus=[0 0 40]/1000;  %  Initial electronic focus

%% %%%%%%%%%%%%%%%%%%%%%%%
%  Define the transducers
Th = xdc_linear_array (N_elements, width, height, kerf, 2, 3, focus);
Th2 = xdc_linear_array (N_elements2, width, height, kerf, 2, 3, focus);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (Th, impulse_response);
xdc_impulse (Th2, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (Th, excitation);

%% %%%%%%%%%%%%%%%%%%%%%%%
%  Define a small phantom with scatterers 设定散射体 200个散射点  随机散射点
N=200; %  Number of scatterers
x_size = 20/1000;  %  Width of phantom [mm]
y_size = 10/1000;  %  Transverse width of phantom [mm]
z_size = 20/1000;  %  Height of phantom [mm]
z_start = 5/1000;  %  Start of phantom surface [mm];

%  Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
positions=[x y z];

%  Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);
%% %%%%%%%%%%%%%%%%%%%%%%%
%  Do the calculation
[v,t]=calc_scat_all (Th, Th2, positions, amp, 1);  %200个散射点，位置随机，散射强度随机amp
%  Plot the individual responses
[N,M]=size(v);
scale=max(max(v));
v=v/scale;
for i=1:M
plot((0:N-1)/fs+t,v(:,i)+i,'b'), hold on
end
hold off





