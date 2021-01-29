% Course: ELEC4700A, LAB 3
% Name: Oritseserundede Eda (Albert)
% CUID: 100993421

% This is a simple 1D electron simulation tool that takes the initial
% position and velocity of an electron at x & v = 0.

% With a force present in time, the position and velocity of the electron
% are updated. 

% The parameters involved in this simulation

e_mass = 9.11E-31; % the mass of an electron (in kg)
e_charge = 1.60217662E-19; % charge of an electron (in coulombs)
e_field = 1E5; % the induced electric field (volt/m)
num_e = 100; % number of electrons 
p_e = 0.05; % probability of electron deflections 
dt = 1E-6; % the timestep used in the simulation 
t_max = 1E-3; % the end time for our simulation 

% The vectors for the postion and velocity of the 100 electrons 

elec_p = zeros(num_e, 1); % electron position 
elec_v = zeros(num_e, 1); % electron velocity 

% Plotting the electron's initial position, velocity and drift velocity
% before the movie displacement of the electron. 
% We would need the plot to be established before going into our loop. 

% the electron initial position 
subplot(3,1,1);
title ('Electron Position'); 
ylabel ('Displacement (m)');
hold on; 

% the electron initial velocity 
subplot(3,1,2);
title ('Electron Velocity');
ylabel ('Velocity (m/s)'); 
hold on; 

% the electron initial drift velocity 
subplot(3,1,3);
title ('Electron Drift Velocity');
ylabel ('Velocity (m/s)');
xlabel ('Time (s)'); 
hold on;

% The simulation loop for modelling the movie of the position and velocity
% of the electron over the range from 0 to 0.001 in steps of 0.000001

for t = 0 : dt : t_max
    
    % Calculate the time steps in terms of the zerod vectors for our
    % position and velocity. 
    elec_p = elec_p + elec_v * dt; 
    elec_v = elec_v + ((e_charge * e_field)/e_mass) * dt;
    
    % Scattering the electrons 
    elec_v(rand(num_e,1) < p_e) = 0;
    
    % The Drift Velocity 
    d_v = mean(elec_v); 
    
    % The electron's changing position, velocity an drift velocity 
    % Position 
    subplot(3,1,1);
    plot(repmat(t,num_e,1), elec_p, 'rX');
    hold on; 
    
    % Velocity 
    subplot(3,1,2);
    plot(repmat(t,num_e,1), elec_v, 'gX'); 
    hold on; 
    
    % Drift velocity
    subplot(3,1,3);
    plot(t,d_v,'bX'); 
    hold on; 
    
    %pause(0.01);
end
