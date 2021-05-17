%% Project 1

% Daniel Bristow
% PHYS 430

%{

This model is NOT corrected for general relativity. This is not a major
issue for simulations of table orbits with very low eccentricity, but this
will be inaccuracte for a chaotic system.

For uses of this program, it is reccomended to try out different values for
the M_J_MULT, M_E_MULT, and M_S_MULT variables. They are multipliers of the
real masses of Jupiter, Earth, and the Sun respectively. Note that the
variabels are redefined in each section.

Also note that due to the small time step for the large amount of time (as
well as the fact that this is some fo my earlier MATLAB code), the program 
is quite memory intensive! Feel free to alter the time step DT, but note
that it must be very small for the simulation to be accurate, else the
bodies may be ejected from orbit in the simulation when in (a classical)
reality they wouldn't.

DT is time step (size of iteration).
KMAX is the number of iterations.

%}

% NOTE: Initial iteration k=1

% Unit of distance: AU
% Unit of mass: kg (not used in calculations, only mass ratios are needed)
% Unit of time: yr

clear;
close all;


%% Example 4.1

% Earth and Sun

% Sun is center of mass

% Earth's distance from the sun is 1 AU
% Earth's linear velocity is 2*pi AU/yr

KMAX = 500;     % iterations
DT = 0.002;     % time step (in years)
X1 = 1;         % inital x
Y1 = 0;         % inital y
VX1 = 0;        % inital x velocity
VY1 = 2*pi;     % inital y velocity

x = zeros(1,KMAX); y = zeros(1,KMAX);

vx = zeros(1,KMAX); vy = zeros(1,KMAX);


x(1,1) = X1; y(1,1) = Y1;

vx(1,1) = VX1; vy(1,1) = VY1;


for k = 1 : KMAX-1
    
    rk = sqrt(x(1,k)^2 + y(1,k)^2); % Distance from sun
    
    vx(1,k+1) = vx(1,k) - ((4*pi^2*x(1,k))/rk^3)*DT;
    vy(1,k+1) = vy(1,k) - ((4*pi^2*y(1,k))/rk^3)*DT;
    
    % Euler-Cromer step
    x(1,k+1) = x(1,k) + vx(1,k+1)*DT;
    y(1,k+1) = y(1,k) + vy(1,k+1)*DT;    
    
end

figure (1)

plot(x,y, '.')

xlim([-1,1])
ylim([-1,1])
pbaspect([1 1 1])

title(strcat('Earth and Sun'))
xlabel('x (AU)')
ylabel('y (AU)')

%% Example 4.2

% Earth, Jupiter, and Sun

% Sun is center of mass

%{

We begin to see a chaotic orbit for Earth when Jupiter is set to 1000 times
its mass. However, because the sun is still the center of mass here, it is 
not an accurate classical simulation.

%}

% Earth's distance from the sun is 1 AU
% Earth's linear velocity is 2*pi AU/yr

clear;

KMAX = 2500000;   % iterations (50000000 takes about a minute!)
DT = 0.000002;    % time step (in years) (no ejection so far at 0.000002)

% EARTH
M_E = 6.0*10^24;
X_E1 = 1;         % inital x
Y_E1 = 0;         % inital y
VX_E1 = 0;        % inital x velocity
VY_E1 = 2*pi;     % inital y velocity

% JUPITER
M_J_MULT = 1000;            % How many times Jupiter mass
M_J = 1.9*10^27*M_J_MULT;
X_J1 = 5.2;                 % inital x (5.2 AU on average)
Y_J1 = 0;                   % inital y
VX_J1 = 0;                  % inital x velocity
VY_J1 = 2*pi*(5.2/12);      % inital y velocity (12 year orbit)

% SUN
M_S = 2.0*10^30;


x_e = zeros(1,KMAX); y_e = zeros(1,KMAX);
vx_e = zeros(1,KMAX); vy_e = zeros(1,KMAX);

x_j = zeros(1,KMAX); y_j = zeros(1,KMAX);
vx_j = zeros(1,KMAX); vy_j = zeros(1,KMAX);


x_e(1,1) = X_E1; y_e(1,1) = Y_E1;
vx_e(1,1) = VX_E1; vy_e(1,1) = VY_E1;

x_j(1,1) = X_J1; y_j(1,1) = Y_J1;
vx_j(1,1) = VX_J1; vy_j(1,1) = VY_J1;


for k = 1 : KMAX-1
    
    rk_e = sqrt(x_e(1,k)^2 + y_e(1,k)^2); % Earth distance from sun
    rk_j = sqrt(x_j(1,k)^2 + y_j(1,k)^2); % Jupiter distance from sun
    rk_e_j = sqrt((x_e(1,k) - x_j(1,k))^2 + (y_e(1,k) - y_j(1,k))^2);
    
    vx_e(1,k+1) = vx_e(1,k) - ((4*pi^2*x_e(1,k))/rk_e^3)*DT - ((4*pi^2*(M_J/M_S)*(x_e(1,k) - x_j(1,k)))/rk_e_j^3)*DT;
    vy_e(1,k+1) = vy_e(1,k) - ((4*pi^2*y_e(1,k))/rk_e^3)*DT - ((4*pi^2*(M_J/M_S)*(y_e(1,k) - y_j(1,k)))/rk_e_j^3)*DT;
    
    vx_j(1,k+1) = vx_j(1,k) - ((4*pi^2*x_j(1,k))/rk_j^3)*DT - ((4*pi^2*(M_E/M_S)*(x_j(1,k) - x_e(1,k)))/rk_e_j^3)*DT;
    vy_j(1,k+1) = vy_j(1,k) - ((4*pi^2*y_j(1,k))/rk_j^3)*DT - ((4*pi^2*(M_E/M_S)*(y_j(1,k) - y_e(1,k)))/rk_e_j^3)*DT;
    
    % Euler-Cromer step
    x_e(1,k+1) = x_e(1,k) + vx_e(1,k+1)*DT;
    y_e(1,k+1) = y_e(1,k) + vy_e(1,k+1)*DT;
    
    x_j(1,k+1) = x_j(1,k) + vx_j(1,k+1)*DT;
    y_j(1,k+1) = y_j(1,k) + vy_j(1,k+1)*DT;    
    
end

figure (2)

plot(x_e,y_e, '.')
hold on
plot(x_j,y_j, '.')
hold off

legend('Earth', 'Jupiter')

xlim([-6,6])
ylim([-6,6])
pbaspect([1 1 1])
grid on

title(strcat('Earth, Jupiter, and Sun (M_{Sun}>>M_{Jupiter}>M_{Earth})'))
xlabel('x (AU)')
ylabel('y (AU)')

%% Exercise 4.16

%{

We now make the consideration for the center of mass of the three body
system being at the origin. The mass of the sun is no longer dominant
enough to be considered the center of mass.

Note that this will still be innacurate being that relativity is not being
considered here. This is still a classical simulation.

To create this simultion, we align all three bodies such that the center of
mass is at the origin. We then set the initial velocities of all three
bodies (including the sun) such that the total momentum of the system is
zero and the center of mass remains at the origin.

%}

clear;

KMAX = 2500000;   % iterations (50000000 takes about a minute!)
DT = 0.00002;    % time step (in years)

% EARTH
M_E_MULT = 1;
M_E = 6.0*10^24*M_E_MULT;
X_E1 = 1;         % inital x
Y_E1 = 0;         % inital y (KEEP AT 0)
VX_E1 = 0;        % inital x velocity (KEEP AT 0)
VY_E1 = 2*pi;     % inital y velocity

% JUPITER
M_J_MULT = 100;
M_J = 1.9*10^27*M_J_MULT;
X_J1 = 5.2;                 % inital x (5.2 AU on average)
Y_J1 = 0;                   % inital y (KEEP AT 0)
VX_J1 = 0;                  % inital x velocity (KEEP AT 0)
VY_J1 = 2*pi*(5.2/12);      % inital y velocity (12 year orbit)

% SUN (position and velocity based off inputs for Earth and Jupiter)
M_S_MULT = 1;
M_S = 2.0*10^30;

% Gravitational constant (AU^3*kg^-1*yr^-2)
G = 4*pi^2/(M_S/M_S_MULT);


x_e = zeros(1,KMAX); y_e = zeros(1,KMAX);
vx_e = zeros(1,KMAX); vy_e = zeros(1,KMAX);

x_j = zeros(1,KMAX); y_j = zeros(1,KMAX);
vx_j = zeros(1,KMAX); vy_j = zeros(1,KMAX);

x_s = zeros(1,KMAX); y_s = zeros(1,KMAX);
vx_s = zeros(1,KMAX); vy_s = zeros(1,KMAX);


x_e(1,1) = X_E1; y_e(1,1) = Y_E1;
vx_e(1,1) = VX_E1; vy_e(1,1) = VY_E1;

x_j(1,1) = X_J1; y_j(1,1) = Y_J1;
vx_j(1,1) = VX_J1; vy_j(1,1) = VY_J1;

x_s(1,1) = -((M_E*X_E1 + M_J*X_J1)/M_S); y_s(1,1) = 0;
vx_s(1,1) = 0; vy_s(1,1) = -((M_E*VY_E1 + M_J*VY_J1)/M_S);


for k = 1 : KMAX-1
    
    rk_e_s = sqrt((x_e(1,k) - x_s(1,k))^2 + (y_e(1,k) - y_s(1,k))^2); % Earth distance from sun
    rk_j_s = sqrt((x_j(1,k) - x_s(1,k))^2 + (y_j(1,k) - y_s(1,k))^2); % Jupiter distance from sun
    rk_e_j = sqrt((x_e(1,k) - x_j(1,k))^2 + (y_e(1,k) - y_j(1,k))^2); % Earth distance from Jupiter
    
    vx_e(1,k+1) = vx_e(1,k) - ((G*M_J*(x_e(1,k) - x_j(1,k)))/rk_e_j^3)*DT - ((G*M_S*(x_e(1,k) - x_s(1,k)))/rk_e_s^3)*DT;
    vy_e(1,k+1) = vy_e(1,k) - ((G*M_J*(y_e(1,k) - y_j(1,k)))/rk_e_j^3)*DT - ((G*M_S*(y_e(1,k) - y_s(1,k)))/rk_e_s^3)*DT;
    
    vx_j(1,k+1) = vx_j(1,k) - ((G*M_E*(x_j(1,k) - x_e(1,k)))/rk_e_j^3)*DT - ((G*M_S*(x_j(1,k) - x_s(1,k)))/rk_j_s^3)*DT;
    vy_j(1,k+1) = vy_j(1,k) - ((G*M_E*(y_j(1,k) - y_e(1,k)))/rk_e_j^3)*DT - ((G*M_S*(y_j(1,k) - y_s(1,k)))/rk_j_s^3)*DT;
    
    vx_s(1,k+1) = vx_s(1,k) - ((G*M_E*(x_s(1,k) - x_e(1,k)))/rk_e_j^3)*DT - ((G*M_J*(x_s(1,k) - x_j(1,k)))/rk_j_s^3)*DT;
    vy_s(1,k+1) = vy_s(1,k) - ((G*M_E*(y_s(1,k) - y_e(1,k)))/rk_e_j^3)*DT - ((G*M_J*(y_s(1,k) - y_j(1,k)))/rk_j_s^3)*DT;
    
    % Euler-Cromer step
    x_e(1,k+1) = x_e(1,k) + vx_e(1,k+1)*DT;
    y_e(1,k+1) = y_e(1,k) + vy_e(1,k+1)*DT;
    
    x_j(1,k+1) = x_j(1,k) + vx_j(1,k+1)*DT;
    y_j(1,k+1) = y_j(1,k) + vy_j(1,k+1)*DT;
    
    x_s(1,k+1) = x_s(1,k) + vx_s(1,k+1)*DT;
    y_s(1,k+1) = y_s(1,k) + vy_s(1,k+1)*DT;  
    
end

figure (3)

plot(x_e,y_e, '.')
hold on
plot(x_j,y_j, '.')
plot(x_s,y_s, '.')
hold off

legend('Earth', 'Jupiter', 'Sun')

xlim([-6,6])
ylim([-6,6])
pbaspect([1 1 1])
grid on

title(strcat('Earth, Jupiter, and Sun'))
xlabel('x (AU)')
ylabel('y (AU)')
