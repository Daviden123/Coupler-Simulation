clear;
clc;
% Define the symbolic variables
syms V1 V2 V3 V4 V5 V6 C_c1 C_c2

% Define constants
f = 13.56e6;  % Frequency in Hz
Vin = 1;      % Input voltage
Z0 = 50;      % Characteristic impedance
L = 500e-9;   % Transmission Line Lumped Element Inductors
C = 220e-12;  % Transmission Line Lumped Element Capacitors
Z_L = 2 * pi * f * L*1j; % Impedance of Lumped Element Inductors
Z_C = -1j*1 /(2 * pi * f * C); % Impedance of Lumped Element Capacitors
X1 = 1 /(2 * pi * f * C_c1);  % Reactance of Coupling Capacitors on the side
X2 = 1 /(2 * pi * f * C_c2);  % Reactance of Coupling Capacitors in the middle

% Define the nodal equations (Nodal Analysis)
eq1 = (V1 - Vin) / Z0 + (V1 - V3) / (-X1*1j) + (V1 - V5) / Z_L == 0;
eq2 = (V5 - V1) / Z_L + (V5) / Z_C + (V5 - V2) / Z_L + (V5 - V6)/ (-X2*1j) == 0;
eq3 = (V2 - V5) / Z_L + (V2 - V4) / (-X1*1j) + (V2) / Z0 == 0;
eq4 = V3 / Z0 + (V3 - V1) / (-X1*1j) + (V3 - V6) / Z_L == 0;
eq5 = (V6 - V3) / Z_L + (V6) / Z_C + (V6 - V4) / Z_L  + (V6 - V5) / (-X2*1j) == 0;
eq6 = (V4 - V6) / Z_L + (V4 - V2) / (-X1*1j) + (V4) / Z0 == 0;

% Collect equations into a system
equations = [eq1; eq2; eq3; eq4; eq5; eq6];

% Define the variables to solve for
variables = [V1 V2 V3 V4 V5 V6];

% Solve the system of equations
solutions = solve(equations, variables);

% Extract the solutions
V1_sol = simplify(solutions.V1);
V2_sol = simplify(solutions.V2);
V3_sol = simplify(solutions.V3);
V4_sol = simplify(solutions.V4);
V5_sol = simplify(solutions.V5);
V6_sol = simplify(solutions.V6);

% Compute the S-parameters
S11 = simplify((V1_sol - Vin) / Vin);
S21 = simplify(2*V2_sol / Vin);
S31 = simplify(2*V3_sol / Vin);
S41 = simplify(2*V4_sol / Vin);

% Compute the magnitudes of the S-parameters
S11_mag = abs(S11);
S21_mag = abs(S21);
S31_mag = abs(S31);
S41_mag = abs(S41);

% Define a list of specific values for available Capacitance in the Market
%C_c_values = [10, 11, 12, 13, 15, 18, 20, 22, 24, 27, 30, 33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91, 100, 110, 120, 130, 150, 180, 200, 220] * 1e-12;
min_C_c_value = 10e-12;
max_C_c_value = 250e-12;
C_c_value_Step = 1e-12;
C_c_values = min_C_c_value:C_c_value_Step:max_C_c_value;

% Initialize lists to store the evaluated magnitudes of the S-parameters
S21_mag_values = [];
S31_mag_values = [];
S41_mag_values = [];

% Loop over the list of C_c values and evaluate the S-parameters
for C_c1_value = C_c_values
    for C_c2_value = C_c_values
        S21_mag_value = double(subs(abs(S21), {C_c1, C_c2}, {C_c1_value, C_c2_value}));
        S31_mag_value = double(subs(abs(S31), {C_c1, C_c2}, {C_c1_value, C_c2_value}));
        S41_mag_value = double(subs(abs(S41), {C_c1, C_c2}, {C_c1_value, C_c2_value}));
        
        S21_mag_values = [S21_mag_values, S21_mag_value];
        S31_mag_values = [S31_mag_values, S31_mag_value];
        S41_mag_values = [S41_mag_values, S41_mag_value];
    end
end
% Convert magnitudes to dB
S21_dB_values = 20 * log10(S21_mag_values);
S31_dB_values = 20 * log10(S31_mag_values);
S41_dB_values = 20 * log10(S41_mag_values);

%% Search S21 dB
Target_S41 = -3; % Desired Coupling Coeff
db_tolerance = 0.1; % Acceptable Coupling Coeff Deviation
Isolation_max = 0; % Init
Output_S41_index = 1; % Init

% Search maximum isolation within tolerance
Search_array = abs(ones(1,length(S41_dB_values))*Target_S41-S41_dB_values);
for i = 1:length(Search_array)
    if((Search_array(i) <= db_tolerance) && (S31_dB_values(i) < Isolation_max))
        Output_S41_index = i;
        Isolation_max = S31_dB_values(i);
    end
end

Opt_CC1 = (C_c_values(fix(Output_S41_index / length(C_c_values))+1))*1e12;
Opt_CC2 = (C_c_values(mod((Output_S41_index-1) , length(C_c_values))+1))*1e12;

% Display the evaluated S-parameters in dB
fprintf('Evaluated S-parameters at Target = %.0f dB, C_c1 = %.0f pF, C_c2 = %.0f pF:\n', Target_S41, Opt_CC1, Opt_CC2);
fprintf('S21 (dB) transmission = %.4f\n', S21_dB_values(Output_S41_index));
fprintf('S31 (dB) isolation = %.4f\n', S31_dB_values(Output_S41_index));
fprintf('S41 (dB) couple = %.4f\n', S41_dB_values(Output_S41_index));
fprintf('\n');