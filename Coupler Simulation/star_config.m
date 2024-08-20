% Define the symbolic variables
syms V1 V2 V3 V4 V5 V6 C_c

% Define constants
f = 13.56e6;  % Frequency in Hz
Vin = 1;      % Input voltage
Z0 = 50;      % Characteristic impedance
L = 500e-9;
C = 220e-12;
Z_L = 2 * pi * f * L*1j;
Z_C = -1j*1 /(2 * pi * f * C);
X = 1 /(2 * pi * f * C_c);  % Reactance

% Define the nodal equations
eq1 = (V1 - Vin) / Z0 + (V1 - V3) / (-X*1j) + (V1 - V5) / Z_L == 0;
eq2 = (V5 - V1) / Z_L + (V5) / Z_C + (V5 - V2) / Z_L == 0;
eq3 = (V2 - V5) / Z_L + (V2 - V4) / (-X*1j) + (V2) / Z0 == 0;
eq4 = V3 / Z0 + (V3 - V1) / (-X*1j) + (V3 - V6) / Z_L == 0;
eq5 = (V6 - V3) / Z_L + (V6) / Z_C + (V6 - V4) / Z_L == 0;
eq6 = (V4 - V6) / Z_L + (V4 - V2) / (-X*1j) + (V4) / Z0 == 0;

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

% Define a list of specific values for C_c
C_c_values = [10, 11, 12, 13, 15, 18, 20, 22, 24, 27, 30, 33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91, 100, 110, 120, 130, 150, 180, 200, 220] * 1e-12;

% Initialize lists to store the evaluated magnitudes of the S-parameters
S21_mag_values = [];
S31_mag_values = [];
S41_mag_values = [];

% Loop over the list of C_c values and evaluate the S-parameters
for C_c_value = C_c_values
    S21_mag_value = double(subs(abs(S21), C_c, C_c_value));
    S31_mag_value = double(subs(abs(S31), C_c, C_c_value));
    S41_mag_value = double(subs(abs(S41), C_c, C_c_value));
    
    S21_mag_values = [S21_mag_values, S21_mag_value];
    S31_mag_values = [S31_mag_values, S31_mag_value];
    S41_mag_values = [S41_mag_values, S41_mag_value];
end

% Convert magnitudes to dB
S21_dB_values = 20 * log10(S21_mag_values);
S31_dB_values = 20 * log10(S31_mag_values);
S41_dB_values = 20 * log10(S41_mag_values);

% Display the evaluated S-parameters in dB
for i = 1:length(C_c_values)
    fprintf('Evaluated S-parameters at C_c = %.2e F:\n', C_c_values(i));
    fprintf('S21 (dB) transmission = %.4f\n', S21_dB_values(i));
    fprintf('S31 (dB) isolation = %.4f\n', S31_dB_values(i));
    fprintf('S41 (dB) couple = %.4f\n', S41_dB_values(i));
    fprintf('\n');
end