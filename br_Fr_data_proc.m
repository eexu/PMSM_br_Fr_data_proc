clc, clear, close all
% optimized the speed 
tic
% Important parameters settings
sim_num = 3;
space_points_num = 2001; %  this should be a odd number
time_points_num = 200; % must match the csv data
space_order_half = floor(space_points_num/2-1); % this number should be less than or equal to 0.5*space_num -1
%space_order_half = 50; % this number should be less than or equal to 0.5*space_num -1
time_order_half = 100; % this number should be less than or equal to half of time_points_num
n_rpm = 3000; % motor speed
time_period = 0.02/10;
plot_on = 0; % 0: not plot 1: plot
p = 10; % pole pairs
N = 400; % Most N br filter. Reconstruction: Most N parameters set 
top_N = 100; % Reconstruciton Top N This num is recommended to be > 10
most_influential_spatial_order = 2;
most_influential_frequency_order = 20;

%-------------------------------------------------------------------------

Fs = 1/(time_period/time_points_num);
ts = time_period/time_points_num;
most_influential_frequency = most_influential_frequency_order*n_rpm/60; % Hz
% Data initilize
Fr_cells = cell(sim_num, 2); % 1st column is amplitude, 2nd is the phase
br_cells = cell(sim_num, 2); % 1st column is amplitude, 2nd is the phase
Fr_most_amp = zeros(sim_num, 1);
Fr_most_pha = zeros(sim_num, 1);
Fr_comp_cells = cell(sim_num, 1);
br_most_N_cells = cell(sim_num, 1);
spatial_pos = space_order_half+1+most_influential_spatial_order;
frequency_pos = round((most_influential_frequency)/(Fs/time_points_num)+1);
one_hundred_pi = 100*pi;

% fft2 calculation of both Fr and br
for i = 1:(2*sim_num)
    file_name = strcat("Calculator Expressions Plot ", num2str(i), ".csv");
    data_csv = csvread(file_name, 1, 0); %#ok<*CSVRD>
    if i <= sim_num % fft2 of Fr
        [Fr_cells{i, 1}, Fr_cells{i, 2}] = fft2_xt(data_csv, time_period, space_points_num, plot_on, space_order_half, time_order_half);
        Fr_most_amp(i) = Fr_cells{i, 1}(spatial_pos, frequency_pos); % (2, 20f1) Fr
        Fr_most_pha(i) = Fr_cells{i, 2}(spatial_pos, frequency_pos); % (2, 20f1) Fr
    else % fft2 of br
        [br_cells{i, 1}, br_cells{i, 2}] = fft2_xt(data_csv, time_period, space_points_num, plot_on, space_order_half, time_order_half);
    end
end

% Fr reconstruction:
% data preparation
t = 0:ts:(time_period-ts);  % time
theta = 0:2*pi/(space_points_num-1):(2*pi-2*pi/(space_points_num-1)); % space
[T, THETA] = meshgrid(t(1:(end)), theta(1:(end)));
nrow = space_points_num - 1;
ncol = time_points_num;
one_over_L = 1/nrow/ncol;
for selected_num = 1:sim_num
%for selected_num = 7:7
    bar7_amp = br_cells{sim_num+selected_num,1};
    bar7_pha = br_cells{sim_num+selected_num,2};
    % bar7_amp(bar7_amp<1e-5) = 0;  % to show the results clearly
    % bar7_pha(bar7_amp<1e-5) = 0;  % to show the results clearly
    vec_sort = sort(bar7_amp(:), 'descend');
    indices = find(ismember(bar7_amp, vec_sort(1:N))); % find the N largest
    [row, col] = find(ismember(bar7_amp, vec_sort(1:N)));
    Large_number = 520e5; % This number is used to add the first element of br

    res = zeros(N, 4);
    res(:, 1) = (row-(space_order_half+1));     % space order
    res(:, 2) = (col-1)*p;                      % frequency order
    res(:, 3) = bar7_amp(indices);              % harmonics amplitude
    res(:, 4) = bar7_pha(indices);              % phase (rad)
    res(:, 5) = res(:, 3).*cos(res(:, 4));      % x value
    res(:, 6) = res(:, 3).*sin(res(:, 4));      % y value
    res_s = sortrows(res, 3, 'descend');

    target_value = zeros(N, 1);
    Fr_most_delta = zeros(N, 2);
    Fr_most_delta(:, 1) = 1:N;

    z = zeros(size(T));
    for j = 1:N
        % if mod(j,50) == 0 % to display the process
        %     disp(j);
        % end
        z = z + res_s(j, 3) * cos(res_s(j, 2)*one_hundred_pi*T - res_s(j, 1)*THETA + res_s(j, 4));
        Z = z.^2/(8*pi*1e-7);   % Fr
        F2 = fft2(Z);
        F2_rad = angle(F2);
        F2 = abs(F2*one_over_L);
        F2_rad_shift = fftshift(F2_rad);
        F2shift = fftshift(F2);
        time_mid = round(ncol/2+1);
        space_mid = round(nrow/2+1);
        % time_order_half  = 50;  % display time order
        % space_order_half = 220; % display space order
        bar3_data_phase = F2_rad_shift((space_mid-space_order_half):(space_mid+space_order_half), (time_mid-time_order_half):(time_mid));
        bar3_data = F2shift((space_mid-space_order_half):(space_mid+space_order_half), (time_mid - time_order_half):(time_mid));

        bar3_data_non_zero_axis = bar3_data;
        bar3_data_non_zero_axis((space_order_half+1), :) = 0;
        bar3_data_non_zero_axis(:,end) = 0;
        bar3_data = bar3_data + bar3_data_non_zero_axis;
        bar3_data = fliplr(bar3_data);   % positive spatial order on the right side
        bar3_data_phase = fliplr(bar3_data_phase);

        % bar3_data(bar3_data<1e-5) = 0;

        target_value(j) = bar3_data(spatial_pos, frequency_pos); 

        if j == 1
            Fr_most_delta(j, 2) = Large_number;
        else
            Fr_most_delta(j, 2) = abs(target_value(j) - target_value(j-1)); % Calculate delta value
        end
    end

    %%
    Fr_most_delta_sorted = sortrows(Fr_most_delta, 2, "descend"); % get sorted data
    top_N_indices = Fr_most_delta_sorted(1:top_N, 1); % get the top N indices

    most_N_influential_comp = [];
    for i = 1:top_N
        most_N_influential_comp = [most_N_influential_comp; res_s(top_N_indices(i), :)];
        vec = find_s_t_order(res_s(:, 1:2), most_influential_spatial_order, most_influential_frequency_order, res_s(top_N_indices(i), 1), res_s(top_N_indices(i), 2));
        most_N_influential_comp = [most_N_influential_comp; res_s(vec, :)];
    end
    most_N_influential_comp = unique(most_N_influential_comp, 'rows');
    br_most_N_cells{selected_num, 1} = sortrows(most_N_influential_comp, 3, "descend");

    % Use the filtered components to re-generate the target Fr and display
    [filtered_row_n, ~] = size(most_N_influential_comp);
    disp(['Selected case: ', num2str(selected_num)]);
    fprintf('------------------------------------------------------------------\n')
    disp(['Simulation Value & Phase of (', num2str(most_influential_spatial_order), ', ', num2str(most_influential_frequency_order), 'f1):'])
    fprintf(2, [num2str(Fr_most_amp(selected_num)), ' Pa \t| ', num2str(Fr_most_pha(selected_num)), ' rad\n']);
    disp(['Reconstruction Value & Phase of (', num2str(most_influential_spatial_order), ', ', num2str(most_influential_frequency_order), 'f1):'])
    % disp('Reconstruction Value & Phase of (-4, 44f1):')
    % Reconstruct data
    most_N_copied = most_N_influential_comp;
    Fr_re_construction_comp = [];
    while(~isempty(most_N_copied))
        tmp_mat = Fr_re_construct(most_N_copied, most_influential_spatial_order, most_influential_frequency_order);
        if(~isempty(tmp_mat))
            Fr_re_construction_comp = [Fr_re_construction_comp; tmp_mat];
        end
        most_N_copied = most_N_copied(2:end, :);
    end
    Fr_re_construction_comp(:, 7) = Fr_re_construction_comp(:, 5).*cos(Fr_re_construction_comp(:, 6));
    Fr_re_construction_comp(:, 8) = Fr_re_construction_comp(:, 5).*sin(Fr_re_construction_comp(:, 6));
    Fr_comp_cells{selected_num, 1} = sortrows(Fr_re_construction_comp, 5, 'descend'); % store data
    sum_x = sum(Fr_re_construction_comp(:, 7));
    sum_y = sum(Fr_re_construction_comp(:, 8));
    amp_re = sqrt(sum_x^2 + sum_y^2)/(8*pi*1e-7);
    pha_re = angle(sum_x + 1j*sum_y);
    fprintf(2, [num2str(amp_re), ' Pa\t| ', num2str(pha_re), ' rad\n']);
    disp('Amplitude error:');
    disp([num2str((amp_re-Fr_most_amp(selected_num))*100/Fr_most_amp(selected_num)), ' %']);
    disp('Phase error:');
    phase_error_re = rad2deg(pha_re-Fr_most_pha(selected_num));
    if abs(phase_error_re) > 180
        if phase_error_re < 0
            phase_error_re = phase_error_re + 360;
        else
            phase_error_re = phase_error_re - 360;
        end
    end
    disp([num2str(phase_error_re), ' deg']);
    disp('Br used number:');
    disp(num2str(size(most_N_influential_comp, 1)));
    disp('Pairs formed number:');
    disp(num2str(size(Fr_re_construction_comp, 1)));
    fprintf('------------------------------------------------------------------\n')

end

% save data
save ana_data.mat Fr_cells br_cells Fr_most_amp Fr_most_pha...
    Fr_comp_cells br_most_N_cells 
toc

function vec = find_s_t_order(data_list, aim_s, aim_t, known_s, known_t)
vec = [];
% case 1
m = aim_s - known_s;
n = aim_t - known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [vec; intersect_mn];
end

% case 2
m = - aim_s + known_s;
n = - aim_t + known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [vec; intersect_mn];
end

% case 3
m = aim_s + known_s;
n = aim_t + known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [vec; intersect_mn];
end
end

function mat = Fr_re_construct(data_list, aim_s, aim_t)
known_s = data_list(1,1);
known_t = data_list(1,2);
mat = [];
% case 1
m = aim_s - known_s;
n = aim_t - known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [data_list(1, 1:2), data_list(intersect_mn, 1:2), data_list(1, 3)*data_list(intersect_mn, 3),...
        data_list(1, 4)+data_list(intersect_mn, 4)];
    mat = [mat; vec];
end

% case 2
m = - aim_s + known_s;
n = - aim_t + known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [data_list(1, 1:2), data_list(intersect_mn, 1:2), data_list(1, 3)*data_list(intersect_mn, 3),...
        data_list(1, 4)-data_list(intersect_mn, 4)];
    mat = [mat; vec];
end

% case 3
m = aim_s + known_s;
n = aim_t + known_t;
indice_m = find(data_list(:,1) == m);
indice_n = find(data_list(:,2) == n);
intersect_mn = intersect(indice_m, indice_n);
if ~isempty(intersect_mn)
    vec = [data_list(1, 1:2), data_list(intersect_mn, 1:2), data_list(1, 3)*data_list(intersect_mn, 3),...
        -data_list(1, 4)+data_list(intersect_mn, 4)];
    mat = [mat; vec];
end

end
