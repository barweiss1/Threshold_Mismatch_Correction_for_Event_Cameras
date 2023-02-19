% %%  Threshold correction algorithm for a single pixel 
%   Asaf Omer and Bar Weiss
%   04.11.22
% 
%
% The idea is to recontruct the signal sampled by the event camera using
% the known threshold and then downsample it with the new desired threshold

%%  Setup

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.


%% define test signal 

% define sampling parameters 
Fs = 1e5; 
Ts =1/Fs;
t = 0:Ts:2;
N = length(t);
orig_thres_pos = 0.23;
orig_thres_neg = -0.19;
f_0 = 1.3;

% define 
%x = sin(2*pi*t*f_0)+2;
%x = sin(2*pi*(t.^2)*f_0)+2;
%x = (1./(t+1)).*sin(2*pi*(t.^2)*f_0);
%x = t.*sin(2*pi*(t.^2)*f_0);
x = sin(2*pi*(t.^2)*f_0);
%x = sawtooth(2*pi*t*f_0);

figure(1);
plot(t,x,'DisplayName','Original Signal','LineWidth',1.5);
title('1D Event Sampling and Reconstruction');
xlabel('time[sec]');
ylabel('amplitude');
legend

[ts,pol] = sample_1D(x,t,orig_thres_pos,orig_thres_neg);

%% filtered square wave

x = square(2*pi*(t)*f_0);

% create filter in Fourier domain 
x_f = fftshift(fft(x));
f = linspace(-Fs/2,Fs/2,N);
fc = 10; % cut-off frequency 
fil = 1./(1+(1j*f/fc));

% plot FFT of signal and filter
figure(4)
subplot(2,1,1);
plot(f,abs(fil));
title('Low Pass Filter FFT');
subplot(2,1,2);
plot(f,abs(x_f));
title('Signal FFT');

% take the inverse FFT to get the time domain filtered signal 
x_fil_f = x_f.*fil; % apply filter in frequency domain
x_fil = ifft(ifftshift(x_fil_f));

% plot filtered signal 
figure(3);
plot(t,x,'DisplayName','Original Signal','LineWidth',1.5);
hold on 
plot(t,x_fil,'DisplayName','Filtered Signal','LineWidth',1.5);
hold off 
title('Square Wave and LPF Filtered Square Wave');
xlabel('time[sec]');
ylabel('amplitude');
legend


% sample 
[ts,pol] = sample_1D(x_fil,t,orig_thres_pos,orig_thres_neg);
% ZOH and FOH reconstruct
x_zoh = ZOH_rec(ts,pol,orig_thres_pos,orig_thres_neg,t);
x_foh = FOH_rec(ts,pol,orig_thres_pos,orig_thres_neg,t);


figure(1);
plot(t,x_fil,'DisplayName','Filtered Signal');
hold on
plot(t,x_zoh,'DisplayName','ZOH Reconstruction');
plot(t,x_foh,'DisplayName','FOH Reconstruction');
hold off
title('1D Event Sampling and Reconstruction');
xlabel('time[sec]');
ylabel('amplitude');
legend

% correction algorithm 
 
new_thres_pos = 0.09;
new_thres_neg = -0.09;

% resample FOH version and ZOH reconstruct for plotting
[resampled_ts,resampled_pol]= sample_1D(x_foh,t,new_thres_pos,new_thres_neg);
resampled_zoh = ZOH_rec(resampled_ts,resampled_pol,new_thres_pos,new_thres_neg,t);

% sample original signal at new threshold 
[true_ts,true_pol]= sample_1D(x_fil,t,new_thres_pos,new_thres_neg);
new_thres_zoh = ZOH_rec(true_ts,true_pol,new_thres_pos,new_thres_neg,t);

figure(2)
plot(t,new_thres_zoh);
hold on
plot(t,resampled_zoh);
hold off
xlabel('time[sec]');
ylabel('amplitude');
legend('Direct Sample','FOH Resample');
grid on



%% ZOH reconstruct

x_zoh = ZOH_rec(ts,pol,orig_thres_pos,orig_thres_neg,t);
figure(1);
hold on
plot(t,x_zoh,'DisplayName','ZOH Reconstruction','LineWidth',1.5);
hold off
legend

%% FOH reconstruct

x_foh = FOH_rec(ts,pol,orig_thres_pos,orig_thres_neg,t);
figure(1);
hold on
plot(t,x_foh,'DisplayName','FOH Reconstruction','LineWidth',1.5);
hold off
legend

%% resampling FOH and comparing to original sampling 

new_thres_pos = 0.09;
new_thres_neg = -0.09;

% resample FOH version and ZOH reconstruct for plotting
[resampled_ts,resampled_pol]= sample_1D(x_foh,t,new_thres_pos,new_thres_neg);
resampled_zoh = ZOH_rec(resampled_ts,resampled_pol,new_thres_pos,new_thres_neg,t);

% sample original signal at new threshold 
[true_ts,true_pol]= sample_1D(x,t,new_thres_pos,new_thres_neg);
new_thres_zoh = ZOH_rec(true_ts,true_pol,new_thres_pos,new_thres_neg,t);

figure(2);
subplot(2,1,1);
plot(t,new_thres_zoh,'LineWidth',1.5);
hold on
plot(t,resampled_zoh,'LineWidth',1.5);
hold off
xlabel('time[sec]');
ylabel('amplitude');
title('Threshold Correction Vs. Ideal Sample');
legend('Ideal Sample','Threshold Correction');
grid on
subplot(2,1,2);
plot(t,x_zoh,'LineWidth',1.5);
xlabel('time[sec]');
ylabel('amplitude');
title('Original Sample');

%% Error Analysis

err_sig = resampled_zoh - new_thres_zoh;

figure(5)
subplot(2,1,1);
plot(t,err_sig);
xlabel('time[sec]');
ylabel('amplitude');
title('Error Signal');
subplot(2,1,2);
stem(true_ts,true_pol);
hold on
stem(resampled_ts,resampled_pol);
hold off
xlabel('time[sec]');
ylabel('Polarity');
title('Timestamps');
legend('Direct Sample','Recontructed');
grid on


MAE = mean(abs(err_sig));
disp("The MAE of the reconstruction is " + MAE);


%% create MAE graph

thres_pos = 0.01:0.0005:0.5;
thres_neg = -thres_pos;
MAE = zeros(1,length(thres_pos));

for i = 1:length(thres_pos)
   
   % sample at the current threshold 
   [ts,pol] = sample_1D(x,t,thres_pos(i),thres_neg(i));
   x_foh = FOH_rec(ts,pol,thres_pos(i),thres_neg(i),t);
   
   % resample FOH version and ZOH reconstruct for plotting
   [resampled_ts,resampled_pol]= sample_1D(x_foh,t,new_thres_pos,new_thres_neg);
   resampled_zoh = ZOH_rec(resampled_ts,resampled_pol,new_thres_pos,new_thres_neg,t);
   
   % calculate error signal 
   err_sig = resampled_zoh - new_thres_zoh;
   MAE(i) = mean(abs(err_sig));
   
end

figure(6)
plot(thres_pos,MAE);
hold on
title('MAE vs. Normalized Original Threshold');
xlabel('$\theta/A$','Interpreter','Latex');
ylabel('$MAE$','Interpreter','Latex');

%% Test Uniformity of the pixels 

% define parameters 
thres_pos = 0.2; % mean positive threshold
thres_neg = -0.2; % mean negative threshold
sigma = 0.02; % std of the threshold mismatch
pixel_x = 10; % number of pixels
pixel_y = 10;

% randomize thresholds of the camera
pos_map = thres_pos + randn(pixel_x,pixel_y)*sigma;
neg_map = thres_neg + randn(pixel_x,pixel_y)*sigma;

% define time difference matrix 
orig_e_sum = zeros(pixel_x,pixel_y);
corrected_e_sum = zeros(pixel_x,pixel_y);

% define cell arrays for ts and pol
ts_cell = cell(pixel_x,pixel_y);
pol_cell = cell(pixel_x,pixel_y);

ts_post_cell= cell(pixel_x,pixel_y);
pol_post_cell = cell(pixel_x,pixel_y);

for i = 1:pixel_x
    for j = 1:pixel_y
        
        % sample signal at given pixel threshold 
        [ts_cell{i,j},pol_cell{i,j}] = sample_1D(x_fil,t,pos_map(i,j),neg_map(i,j));
        
        % correct with FOH algorithm 
        x_foh = FOH_rec(ts_cell{i,j},pol_cell{i,j},pos_map(i,j),neg_map(i,j),t);
        [ts_post_cell{i,j},pol_post_cell{i,j}]= sample_1D(x_foh,t,thres_pos,thres_neg);
        
        % plot reconstruction
        resampled_zoh = ZOH_rec(ts_post_cell{i,j},pol_post_cell{i,j},thres_pos,thres_neg,t);
        orig_zoh = ZOH_rec(ts_cell{i,j},pol_cell{i,j},pos_map(i,j),neg_map(i,j),t);
        figure(5);
        subplot(2,1,1);
        plot(t,resampled_zoh);
        hold on
        title('Post Correction ZOH signals');
        xlabel('time[sec]');
        ylabel('amplitude');
        subplot(2,1,2);
        plot(t,orig_zoh);
        hold on
        xlabel('time[sec]');
        ylabel('amplitude');
        title('Pre Correction ZOH signals');
        
        
        % calculate the mean of the difference between timestamps
        orig_e_sum(i,j) = numel(pol);
        corrected_e_sum(i,j) = numel(resampled_pol);
        
    end
end

%% first and last timestamp analysis 

% define times of the start and end of change 
Tstart = 1.4; % s
Tend = 1.8; % s

[fst_sigma_pre,fst_sigma_post,lst_sigma_pre,lst_sigma_post] = timestamp_error(ts_cell,ts_post_cell,Tstart,Tend);
disp("The Standard Deviation of the first timestamp is " + num2str(fst_sigma_pre,'%e') + " before the correction and " + num2str(fst_sigma_post,'%e') +" after the correction");
disp("The Standard Deviation of the last timestamp is " + num2str(lst_sigma_pre,'%e') + " before the correction and " + num2str(lst_sigma_post,'%e') +" after the correction");

%% plot time difference 

orig_std = std(orig_e_sum,1,'all')/mean(orig_e_sum,'all');
corrected_std = std(corrected_e_sum,1,'all')/mean(corrected_e_sum,'all');


figure(7)
subplot(1,2,1);
plot(orig_e_sum)
colorbar
title('event number before correction')
subplot(1,2,2);
plot(corrected_e_sum)
colorbar
title('event number after correction')


%% Test Performance with 1% estimation error

% define parameters 
thres_pos = 0.2; % mean positive threshold
thres_neg = -0.2; % mean negative threshold
sigma = 0.02; % std of the threshold mismatch
pixel_x = 10; % number of pixels
pixel_y = 10;

% randomize thresholds of the camera
% pos_map = thres_pos + randn(pixel_x,pixel_y)*sigma;
% neg_map = thres_neg + randn(pixel_x,pixel_y)*sigma;

% shift estimated thresholds with 5% error 
pos_map_est = pos_map + randn(pixel_x,pixel_y)*sigma/10;
neg_map_est = neg_map + randn(pixel_x,pixel_y)*sigma/10;

% define time difference matrix 
orig_e_sum = zeros(pixel_x,pixel_y);
corrected_e_sum = zeros(pixel_x,pixel_y);

% define cell arrays for ts and pol
ts_cell = cell(pixel_x,pixel_y);
pol_cell = cell(pixel_x,pixel_y);

ts_post_cell= cell(pixel_x,pixel_y);
pol_post_cell = cell(pixel_x,pixel_y);

for i = 1:pixel_x
    for j = 1:pixel_y
        
        % sample signal at given pixel threshold 
        [ts_cell{i,j},pol_cell{i,j}] = sample_1D(x_fil,t,pos_map(i,j),neg_map(i,j));
        
        % correct with FOH algorithm based on estimated thresholds
        x_foh = FOH_rec(ts_cell{i,j},pol_cell{i,j},pos_map_est(i,j),neg_map_est(i,j),t);
        [ts_post_cell{i,j},pol_post_cell{i,j}]= sample_1D(x_foh,t,thres_pos,thres_neg);
        
        % plot reconstruction
        resampled_zoh = ZOH_rec(ts_post_cell{i,j},pol_post_cell{i,j},thres_pos,thres_neg,t);
        orig_zoh = ZOH_rec(ts_cell{i,j},pol_cell{i,j},pos_map(i,j),neg_map(i,j),t);
        figure(6);
        subplot(2,1,1);
        plot(t,resampled_zoh);
        hold on
        xlabel('time[sec]');
        ylabel('amplitude');
        title('Post Correction ZOH signals');
        subplot(2,1,2);
        plot(t,orig_zoh);
        hold on
        xlabel('time[sec]');
        ylabel('amplitude');
        title('Pre Correction ZOH signals');

        
        % calculate the mean of the difference between timestamps
        orig_e_sum(i,j) = numel(pol);
        corrected_e_sum(i,j) = numel(resampled_pol);
        
    end
end

%% first and last timestamp analysis 

% define times of the start and end of change 
Tstart = 1.4; % s
Tend = 1.8; % s

[fst_sigma_pre,fst_sigma_post,lst_sigma_pre,lst_sigma_post] = timestamp_error(ts_cell,ts_post_cell,Tstart,Tend);
disp("The Standard Deviation of the first timestamp is " + num2str(fst_sigma_pre,'%e') + " before the correction and " + num2str(fst_sigma_post,'%e') +" after the correction");
disp("The Standard Deviation of the last timestamp is " + num2str(lst_sigma_pre,'%e') + " before the correction and " + num2str(lst_sigma_post,'%e') +" after the correction");

