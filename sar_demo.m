% Dependencies targetgenerator3.m, Image processing toolbox, 
clear;                                                                                   
close all;                                                                               
       
c = 3e8;
% SAR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = 200; %meters
x0 = 20e3;
r0 = sqrt(x0^2 + y0^2); % (m)                                                                         
lambda = 0.03; % (m)                                                                     
w = 50; %(m)width of region to be imaged                                                 
l = 50; %(m)length of region to be imaged                                                
V = 50; %(m/s) velocity of platform   
BW = 300e6; %LFM Bandwidth
tau_p = 50e-6; %uncompressed pulse width
tau_res = 1/BW;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_lfm = BW / tau_p;
                                                                                         
del_y = 1; % (m) desired resolution                                                      
theta_B = del_y / r0; %beam width of SAR array                                           
L = lambda / theta_B; % Required length of SAR array.                                    
T = 50e-3; % Need T < ( lambda * r0) / (2 * V * w) to prevent aliasing                   
% Using the parameters above, T < 120ms                                                                                                                                     
T_L = L/V; %How long the SAR collects data                                             
K_L = T_L/(2*T);                                                                     
% We need a minimum of 2 * K_L + 1 samples, we also want a power of 2                                                                                                         
K_L = pow2( ceil(log(2*K_L+1)/log(2))); %increase K to nearest power of 2                                                                                                                              
w_actual = (lambda * r0) / (2 * V * T); %Actual width of imaged area, need to trim output
% and get rid of anything outside of +-(w/2) (+-25m in this case)                        
L_actual = (K_L - 1) * V * T;                                                            
del_y_actual = (r0 * lambda) / (2 * L_actual); %actual resolution of SAR image                                                                                            
                                                                                                                         
a0 = (2 * V^2) / (lambda * r0);                                                          
k = -K_L/2 : K_L/2 - 1;
r_res = tau_res * c / 2;

% Popup Dialog Window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = listdlg("ListString", {'9 Hard coded scatters', 'Circle', 'Spiral', 'Large Array', 'Load png image', 'HW12 Scatterers'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(choice) > 1 | length(choice) == 0
    choice = 1;
end

switch choice % Hard Coded Scatterers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        Ns = 9; %number of simulated scatterers    
        Ps = [1, 1, 1, 1, 1, 1, 1, 1, 1];
        x = [0, 0, 0, 20, 20, 20, -20, -20, -20];
        y = [0, -20, 20, 0, -20, 20, 0, -20, 20];
        titlestr = '9 Scatters, (±20m, ±20m)';

    case 2 % Circle Scatterers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ns = 500;
        Ps = ones(1,Ns);
        x = (0.5*l) * ones(1,Ns) .* sin(1:Ns);
        y = (0.5*w) * ones(1,Ns) .* cos(1:Ns);
        titlestr = "500 Scatters, Circle with radius: " + num2str(0.5*w);

    case 3 % Spiral Scatterers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ns = 500;
        spirals = 60;
        theta = 0:1/Ns:4*pi;
        phi = 0:1/Ns:4*pi;
        radius_x = 0.75 * l;
        radius_y = 0.75 * w;
        Ps = ones(1,Ns);
        x = radius_x * sin(theta) .* cos(spirals * theta);
        y = radius_y * sin(theta) .* sin(spirals * theta);
        titlestr = '500 Scatters, Spiral';
        
    case 4 % Array of Scatterers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ns = 400;
        to = linspace(1, sqrt(Ns), sqrt(Ns));
        ti = linspace(1, sqrt(Ns), Ns);
        Ps = ones(1,Ns);
        x = ((w/sqrt(Ns)) * (repmat(1:sqrt(Ns), [1 sqrt(Ns)])) - w/2);
        y = ((l/sqrt(Ns)) * (interp1(to, 1:sqrt(Ns), ti)) - l/2);
        titlestr = "400 Scatters, square array, spacing:" + num2str(w/sqrt(Ns))+" (m)";

    case 5 % Load image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_scatters_cnt = floor(w/(2*r_res)); %spacing scatters out at twice range resolution of radar
        y_scatters_cnt = floor(l/(2*r_res));
        Ns = x_scatters_cnt * y_scatters_cnt;
        [file, path] = uigetfile('*.png;*.jpg;*.jpeg;*.tif');
        photo = imread([path file]);
        photo_gray = (rgb2gray(photo));
        photo_rotated = rot90(photo_gray, 3);
        photo_scaled = imresize(photo_rotated,[x_scatters_cnt y_scatters_cnt]);

        to = linspace(1, sqrt(Ns), sqrt(Ns));
        ti = linspace(1, sqrt(Ns), Ns);

        Ps = double(reshape(photo_scaled', [1, y_scatters_cnt*x_scatters_cnt]));
        x = repmat(1:x_scatters_cnt, [1 x_scatters_cnt]) - w/2;
        y = interp1(to, 1:y_scatters_cnt, ti) - l/2;
        titlestr = num2str(floor(Ns)) + " Scatters, image converted to scatters";

    case 6 % Scatters from HW12
        Ns = 3; %number of simulated scatterers    
        Ps = [1, 1, 1];
        x = [0, 2, -2];
        y = [0, -2, 2];
        titlestr = '3 Scatters, (0,0), (2,-2), (-2,2)';        
   otherwise
        Ns = 1; %number of simulated scatterers    
        Ps = [1];
        x = [0];
        y = [0];
        titlestr = 'Single scatter, center of the scene';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdg = (2*V*y0)/(lambda*r0); %gross doppler removal, see slide 52

rmax = sqrt((x0+l/2)^2 + (abs(y0)+w/2+L_actual/2)^2);
tau_min = 2*(x0 - l/2)/c;
tau_max = 2*rmax/c;
PRI = 50e-3;

tau_r0 = 2*r0/c;

m_min = (tau_min-tau_r0) / tau_res;
m_max = (tau_max-tau_r0) / tau_res;
m = m_min:m_max;

% Form the baseband signal from the scatterers defined above
V_BB = zeros(length(m), length(k));
tau = repmat((tau_r0 + m * tau_res)', [1 length(k)]);

% Form voltages/samples out of matched filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Ns %loop over each scatterer, 
    t = T*k + tau_r0 + m'*tau_res; %operations in k and m dims are vectorized
    ri = sqrt( (x0+x(i))^2 + (y0+y(i)-V*t).^2 );
    bb = targetGenerator3(Ps(i), fdg, t, ri, lambda, tau_p, tau, alpha_lfm);
    V_BB = V_BB + bb;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V_BB = V_BB / max(max(V_BB)); %Normalize V_BB- Not required when using imagesc

% Perform RCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_p_cell = 1/BW;
f = k' * 1/(K_L*tau_p_cell); %fft tap frequencies
delta_tau = 2*(sqrt(x0^2 + (y0-V*k*T).^2) - x0)/c;
V_BB_fft = fftshift(fft(V_BB, K_L),1);
RCMC = exp(1j*2*pi*f*delta_tau);
V_BB_corrected = ifft(fftshift(RCMC.*V_BB_fft,1));
V_BB_corrected = V_BB_corrected(1:length(m),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image Sharpening Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v_h = exp(1j*2*pi*(V^2/(lambda*r0))*(k*T).^2);
% v_h = repmat(v_h, [length(m) 1]); %attempted fix for v_h, does not work
v_deltaf = exp(1j*2*pi*( (2*y0*V)/(lambda*r0) * (m'*(c*tau_res/2)/r0)) * k*T);
v_q = exp(-1j*pi*(2*V^2 / (lambda*r0) * m'*(c*tau_res/2)/r0) * (k*T).^2);
%V_BB_corrected = V_BB_corrected .* v_q .* v_deltaf; %UNCOMMENT FOR IMAGE SHARPENING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sar_image = zeros(length(m), length(k));
for i=1:length(m)
    quadratic_phase_correction = (exp(1j * pi * a0 * T^2 * k.^2));%calculate quadratic phase correction
    V_o = fftshift(fft(V_BB_corrected(i,:) .* quadratic_phase_correction)); %take fft across m for each k
    sar_image(i,:) = V_o;
end

sar_image = abs(sar_image);
%sar_image = abs(max(max(sar_image)) - sar_image) / abs(max(max(sar_image)));

del_f = 1 / (T * K_L); % (Hz) frequency resolution 
f = (-K_L/2 : (K_L/2 - 1) ) * del_f;
x = m * c * tau_res / 2;
y = ((lambda*r0*f)/(2*V));

delta_y_corr = (y0 * m * (c*tau_res/2)) / r0;
%y = (2*V/(lambda*r0))*delta_y_corr;

figure
imagesc(x,y,sar_image);
colormap gray;
set(gca, 'Ydir', 'normal');
axis off;

% title(titlestr);
% grid on;
% xlim([-w/2 w/2]);
% ylim([-l/2 l/2]);

% figure
% blurred_image = imgaussfilt(sar_image);
% imagesc(x,y,blurred_image);
% colormap gray;
% set(gca, 'Ydir', 'normal');

% figure
% montage({sar_image, blurred_image});



