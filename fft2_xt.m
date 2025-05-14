function [bar3_data, bar3_data_phase] = fft2_xt(data_raw, period_t, space_points_num, plot_on,space_order_half, time_order_half)
% INPUT:
% data_raw:         read from csv file of Maxwell output
% period_t:         period of the analysis time
% space_points_num: space points number set in Maxwell 3D plot
% plot_on:          to plot or not. 1: plot, 0: not plot
% version:          1.1
% bug fixed:        exceed index problem
% date:             2024.11.20

time_points_num = length(data_raw)/space_points_num;
ts = period_t/(time_points_num-1);
Fs = 1/ts;

Fr_mat = reshape(data_raw(:,3), space_points_num, time_points_num);
[nr, nc] = size(Fr_mat);
% garantee the even row number and column number
if mod(nr, 2) == 1
    Fr_mat = Fr_mat(1:(end-1),:);
end
if mod(nc, 2) == 1
    Fr_mat = Fr_mat(:, 1:(end-1));
end

F2 = fft2(Fr_mat);
[nrow, ncol] = size(F2);
L = nrow * ncol;
F2_rad = angle(F2);
F2 = abs(F2/L);
F2shift_rad = fftshift(F2_rad);
F2shift = fftshift(F2);
time_mid = round(ncol/2+1);
space_mid = round(nrow/2+1);

%time_order_half  = 30;
%space_order_half = 20;
bar3_data = F2shift((space_mid-space_order_half):(space_mid+space_order_half), (time_mid-time_order_half):(time_mid));
bar3_data_phase = F2shift_rad((space_mid-space_order_half):(space_mid+space_order_half), (time_mid-time_order_half):(time_mid));
bar3_data_non_zero_axis = bar3_data;
bar3_data_non_zero_axis((space_order_half+1),:) = 0;
bar3_data_non_zero_axis(:,end) = 0;
bar3_data = bar3_data + bar3_data_non_zero_axis;
bar3_data = fliplr(bar3_data); % positive spatial order on the right side
bar3_data_phase = fliplr(bar3_data_phase);

% plot
if plot_on
    f = (Fs/ncol)*(0:5:(time_order_half));
    space_order = -space_order_half:5:space_order_half+1;
    b3 = bar3(bar3_data);
    set(gca, 'XTick', 1:5:(time_order_half+1), 'XTickLabel', f);
    set(gca, 'YTick',1:5:space_order_half*2+1 ,'YTickLabel', space_order);
    % bar3 feature:
    % For an m-by-n matrix, bar3 plots the bars on an x-axis ranging from 1 to n and a y-axis ranging from 1 to m.
    xlabel('Frequency (Hz)');
    ylabel('Spatial order');
end

end
