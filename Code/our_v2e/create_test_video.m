% 
% The code below generates sythetic videos at several patterns to be
% inputted to the v2e tool and then used for threhshold estimation.
%
% 


clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.
fontSize = 14;
   
% Change the current folder to the folder of this m-file.
% (The line of code below is from Brett Shoelson of The Mathworks.)
workingDir = fileparts(matlab.desktop.editor.getActiveFilename);
cd(workingDir);

%% Screen lighting up and down
M = 240;
N = 180;
cycles = 10;
len_sec = 1;
T = 25*len_sec;
lin_pattern = (1:T)/T;
vid = zeros(N,M,2*T*cycles);

frames_no = 2400;
vid = zeros(N,M,frames_no);
frame_rate = 1200;

for i = 1:2:cycles
    for t =1:T
       vid(:,:,(i-1)*T + t) = lin_pattern(t);
       vid(:,:,((i)*T) + t) = lin_pattern(T-t+1);
    end
end

outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/flicker_video.avi'));
outputCustomVideo.FrameRate = 20;
open(outputCustomVideo)

for t = 1:T*cycles
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/flicker_video.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);


%% Moving vertical bars through screen
M = 240;
N = 180;

bar_set_width = 10;
bars_number = M/bar_set_width;

all_bars = [];
bar_set = [zeros(N,bar_set_width/2) ones(N,bar_set_width/2)];

for i = 1:bars_number
    all_bars = [all_bars bar_set];
end

frames_no = 500;
vid = zeros(N,M,frames_no);

outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/bars_vertical_video.avi'));
outputCustomVideo.FrameRate = 30;
open(outputCustomVideo)

for t = 1:frames_no
    vid(:,:,t) = all_bars;
    all_bars = circshift(all_bars,1,2);
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/bars_vertical_video.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);

%% Moving vertical bars through screen
M = 240;
N = 180;

bar_set_width = 10;
bars_number = N/bar_set_width;

all_bars = [];
bar_set = [zeros(bar_set_width/2 , M); ones(bar_set_width/2, M)];

for i = 1:bars_number
    all_bars = [all_bars; bar_set];
end

frames_no = 2400;
vid = zeros(N,M,frames_no);

outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/bars_horizontal_video.avi'));
outputCustomVideo.FrameRate = 800;
open(outputCustomVideo)

for t = 1:frames_no
    vid(:,:,t) = all_bars;
    all_bars = circshift(all_bars,1,1);
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/bars_horizontal_video.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);

%%

%% Moving Sine

M = 240;
N = 180;
y = ones(N, 1);
x = linspace(1,M,M);
full_sine = y.*sin(0.1*x);

frames_no = 1000;
vid = zeros(N,M,frames_no);
frame_rate = 120;
outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/moving_sine_video.avi'));
outputCustomVideo.FrameRate = frame_rate ;
open(outputCustomVideo)

for t = 1:frames_no
    vid(:,:,t) = full_sine;
    full_sine = y.*sin(0.1*x-(5*t/frame_rate));
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/moving_sine_video.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);

%% Moving linear Gradient

M = 240;
N = 180;
y = ones(N, 1);
x = linspace(1,M,M);
wavelen = 100;
velocity = 13;
x_t = x/wavelen;
lin_grad = y.*(x_t-floor(x_t));

frames_no = 8000;
vid = zeros(N,M,frames_no);
frame_rate = 1000;
outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/linear_gradient_video_8s.avi'));
outputCustomVideo.FrameRate = frame_rate ;
open(outputCustomVideo)

for t = 1:frames_no
    vid(:,:,t) = lin_grad;
    x_t = x/wavelen - velocity*t/frame_rate;
    sign_x = (-1).^(floor(x_t));
    x_grad = (1-sign_x)/2+(x_t-floor(x_t)).*sign_x;
    lin_grad = y.*x_grad;
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/linear_gradient_video_8s.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);


%% Flicker

M = 240;
N = 180;
velocity = 13;
flicker = zeros(N,M);

frames_no = 8000;
vid = zeros(N,M,frames_no);
frame_rate = 1000;
outputCustomVideo = VideoWriter(fullfile(workingDir,'generated videos/flicker_video.avi'));
outputCustomVideo.FrameRate = frame_rate ;
open(outputCustomVideo)

for t = 1:frames_no
    vid(:,:,t) = flicker;
    x = velocity*t/frame_rate;
    sign_x = (-1).^(floor(x));
    flicker = (1-sign_x)/2+(x-floor(x)).*sign_x*ones(N,M);
    writeVideo(outputCustomVideo,mat2gray(vid(:,:,t)));
end

close(outputCustomVideo)

CustomMovieFullFileName = fullfile(workingDir,'generated videos/flicker_video.avi');
custom_video = VideoReader(CustomMovieFullFileName,'CurrentTime',0);
