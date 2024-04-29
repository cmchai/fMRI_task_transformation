% Simulating ramping-up activity within a finite interval
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 22nd, April, 2024

close all; clear all
% importing hrf
TR = 1.78;  % in second
[hrf, p] = spm_hrf(TR);

time_hrf = (0:length(hrf)- 1)*TR;

figure;
plot(time_hrf, hrf, 'b-o', "LineWidth", 1.5)
yline(0, '--')

% Define the time interval from CTI onset
interval_start = 0 + 0.5;     % seconds
interval_end = 8.75;    % the langest interval in seconds

% Define the time vector within the interval
time_points = interval_start:TR:interval_end;

% Calculate the slope and intercept for a linearly increasing trend
slope = 1; % or 1 / (interval_end - interval_start);
slope2 = 0.4;
intercept = interval_start;

% Generate the linear trend time course
linear_trend = slope * (time_points - interval_start) + intercept;
linear_trend2 = slope2 * (time_points - interval_start) + intercept;
qua_trend = linear_trend;
qua_trend(4) = qua_trend(2);
qua_trend(5) = qua_trend(1);

% Convolve linear trend with HRF
convolved_trend = conv(linear_trend, hrf);
convolved_trend2 = conv(linear_trend2, hrf);
convolved_qua = conv(qua_trend, hrf);
time_conv = (0:length(convolved_trend)- 1)*TR;

% cropped convolved linear trend
crop_convolved_trend = convolved_trend(1:length(time_points));
crop_convolved_trend2 = convolved_trend2(1:length(time_points));
crop_convolved_qua = convolved_qua(1:length(time_points));


figure;
subplot(2, 2, 1);
plot(time_points, linear_trend, 'b-o', 'LineWidth', 2);
hold on
plot(time_points, linear_trend2, 'm-o', 'LineWidth', 2);
hold on
plot(time_points, qua_trend, 'r-o', 'LineWidth', 2);
xline(5,'--');
xline(interval_end,'--');
xlabel('Time (seconds)');
ylabel('neural activity');
title('Hypothesized neural activity');

subplot(2, 2, 2);
plot(time_hrf, hrf, 'k-o', 'LineWidth', 2);
xline(5,'--');
xline(interval_end,'--');
yline(0, '--');
xlabel('Time (seconds)');
ylabel('BOLD signal');
title('Canonical HRF');

subplot(2, 2, 3);
plot(time_points, crop_convolved_trend, 'b-o', 'LineWidth', 2);
hold on
plot(time_points, crop_convolved_trend2, 'm-o', 'LineWidth', 2);
hold on
plot(time_points, crop_convolved_qua, 'r-o', 'LineWidth', 2);
xline(5,'--');
xline(interval_end,'--');
xlabel('Time (seconds)');
ylabel('BOLD signal');
title('BOLD within CTI after convolving with HRF');

subplot(2, 2, 4);
plot(time_conv, convolved_trend, 'b-o', 'LineWidth', 2);
hold on
plot(time_conv, convolved_trend2, 'm-o', 'LineWidth', 2);
hold on
plot(time_conv, convolved_qua, 'r-o', 'LineWidth', 2);
xline(5,'--');
xline(interval_end,'--');
yline(0, '--');
xlabel('Time (seconds)');
ylabel('BOLD signal');
title('complete BOLD response after convolving with HRF');