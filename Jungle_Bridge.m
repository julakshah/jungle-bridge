%% Loading the data
fpath = "./RubberBandTemplate.csv";
my_table = readtable(fpath);
disp(my_table)

row_range = 1:12; % all the rows we need.
col_range = 3:6; % all the columns we need.
disp(my_table(row_range, col_range))
data_mat = table2array(my_table(row_range, col_range));
disp(data_mat)

%% Exercise 22.1

gravity = 9.81;
count_rubberband = 1; % Selecting rubber band #1

% finding the indices for the row/col values.
row_mass = (count_rubberband * 2) - 1; % DONT CHANGE THIS.
row_length = 2; % REPLACE THIS STUFF WITH A FORMULA OF SOME KIND (BEN TARR REFERENCE)

% collecting all of the mass + length values.
values_mass = data_mat(row_mass,:) / 1000; % divided by 1000 for conversion to kg.
values_force = gravity * values_mass; % force acting downward on each mass.

values_length = data_mat(row_length,:) / 100; % divided by 100 for conversion to meters.
transposed_length = values_length'; % Just the lengths in column form.

% constructing the transposed length matrix.
A = [transposed_length, ones(size(transposed_length))];
disp(A);

%% Exercise 22.2

transposed_force = values_force';
LOBF_RB1 = (A' * A) \ (A' * transposed_force);
disp(LOBF_RB1);

k = LOBF_RB1(1); % This is also slope (m).
b = LOBF_RB1(2); % this is just y-int.

natural_length = ((-b) / k);

% Getting the line of best fit.
coefficients = polyfit(values_mass, values_length, 1);
y_fit = polyval(coefficients, values_mass);

% Plotting the points against the line of best fit.
figure;
plot(values_mass, values_length, '*'); % the points.
hold on; grid on;
plot(values_mass, y_fit); % best-fit line.
xlabel('Mass (kilograms)');
ylabel('Length / Displacement (Meters)');
title('Deformed Length Against Mass- Rubber Band 1')
hold off;
 

%% Exercise 22.3

% creating the empty arrays to store our outputs in.
rubber_band_labels = strings(6,1);
k_values = zeros(6,1);
L0_values = zeros(6,1);

for i = 1:6
    row_mass = 2*i - 1;
    row_len = 2*i;
    
    % accounting for mass, force of gravity, and lengths.
    mass_kg = data_mat(row_mass, :) / 1000;
    force_gravity = gravity * mass_kg;
    len_m = data_mat(row_len, :) / 100;
    
    % line of best fit calculations.
    coeffs = polyfit(mass_kg, len_m, 1);
    y_fit = polyval(coeffs, mass_kg);

    % % plotting the actual line of best fit with every loop.
    % figure;
    % plot(mass_kg, len_m, '*'); hold on; grid on;
    % plot(mass_kg, y_fit);
    % xlabel('Mass (kg)');
    % ylabel('Length (m)');
    % title(['Rubber Band ', num2str(i), ' — Length vs. Mass']);

    % Store values
    rubber_band_labels(i) = ['Rubber Band ', num2str(i)];
    k_values(i) = coeffs(1);
    L0_values(i) = -coeffs(2) / coeffs(1);
end

% visualizing the final table, contains spring constant and natural length
% -> MAKE SURE TO ASK ABOUT THIS TMRW IN QEA, values being negative derive
% from the -b/m formula.
results_table = table(rubber_band_labels, k_values, L0_values, ...
    'VariableNames', {'Rubber Bands', 'Spring Constant (k)', 'Natural Length (l0)'});
disp(results_table);

%% 
% import the data
fpath = "./RubberBandTemplate.csv";
my_table = readtable(fpath);
disp(my_table)


% initialize the system parameters
% which contains parameters describing behavior/measurements of bridge
% param_struct.r0 = [x_0;y_0]: coordinates of leftmost vertex
% param_struct.rn = [x_n;y_n]: coordinates of rightmost vertex
% param_struct.num_links: number of rubber bands in bridge
% param_struct.k_list = [k_1;...;k_n]: list of stiffnesses
% param_struct.l0_list = [l0_1;...;l0_n]: list of natural lengths
% param_struct.m_list = [m_1;...;m_(n-1)]: list of weight masses
% param_struct.g = 9.8 m/sec^2: gravitational acceleration
param_struct = struct();
param_struct.r0 = %YOUR CODE HERE
param_struct.rn = %YOUR CODE HERE
param_struct.num_links = %YOUR CODE HERE
param_struct.k_list = %stiffnesses (N/m), %YOUR CODE HERE
param_struct.l0_list = %natural lengths (meters), %YOUR CODE HERE
param_struct.m_list = %mass list (kg), %YOUR CODE HERE
param_struct.g = %gravitational acceleration (m/s^2) %YOUR CODE HERE

