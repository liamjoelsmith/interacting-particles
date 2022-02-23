%%% Basic simulation of trajectories of a single particle in a 2-D box
hold on

T_i = 300; % Temperature in Kelvin
k = 1.38e-23;  % Boltzmann constant
m = 1e-26;  % mass of a particle in kg
radius = 0.01; % particle radius
L = radius*100; % length of cube side (meters)

% randomly position the particle
x_i = rand; 
y_i = rand;

% select initial particle vecocity components
vx_i = 400;   % initial x vecocity in m/s
vy_i = 200;   % initial y vecocity in m/s

Tau = L/min(vx_i,vy_i); % time for particle to traverse box
tmax = Tau*2.5; % total simulation time
tres = Tau/500; % time resolution
N_time = round(tmax/tres); % number of time steps in simulation

% set initial velocity of particle
vx = vx_i;
vy = vy_i;

% initialise particle position/s
x = zeros(1,N_time); 
x(:,1) = x_i;
y = zeros(1,N_time);
y(:,1) = y_i;

% main simulation loop
for i = 2:N_time 
    % update position
    x(i) = x(i-1) + vx*tres;
    y(i) = y(i-1) + vy*tres;
    
    % check if outside box, and reverse direction if so
    ind = find(x(:,i)>L);  
    vx(ind) = -vx(ind);       
    
    ind = find(x(:,i)<0);
    vx(ind) = -vx(ind);
    
    ind = find(y(:,i)>L);
    vy(ind) = -vy(ind);
    
    ind = find(y(:,i)<0);
    vy(ind) = -vy(ind);
end 

% PLOTTING
t = linspace(0,2*pi);  % full circle revolution, used for 'fill' below
% plot the position of the particle for t>=0
for i = 2:N_time-1
    if round(i/5) == i/5  % only plot every 5th frame, to speed up simulation
        cla
        axis equal
        fill(x(1,i)+radius*cos(t), y(1,i)+radius*sin(t), [1 0 0]);
        axis([0 L 0 L])
        f = getframe;
        im(:,:,1,i/5) = rgb2ind(f.cdata,map,'nodither'); % save the time steps of the plot in a movie
    end
end

imwrite(im,map,'particle_trajectory.gif','DelayTime',0,'LoopCount',inf) % save the movie as a gif
