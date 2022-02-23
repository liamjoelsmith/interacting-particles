%%% Basic simulation of trajectories of a single particle in a 2-D box %%%
clear all
clc
set(gca,'TickLabelInterpreter','latex')

% Constants
T_i = 300;  % Temperature in Kelvin
k = 1.38e-23;  % Boltzmann constant
m = 1e-26;  % Mass of a particle in kg
% Lennard-Jones values for Argon
epsilon = 1.656*10^(-21);
sigma = 0.34*10^(-9);
radius = sigma*2^(1/6);  % particle approx radius (argon)
L = sigma*100;  % length of cube side (meters)
A = 4*epsilon*((sigma)^(12));
B = 4*epsilon*((sigma)^6);

p = 200;  % number of particles
sims = 1;  % number of simulations

% initially space particles evenly in a square-like lattice
N_div = ceil(sqrt(p));
x_0 = transpose(repmat(1:N_div,[1 N_div]));
y_0 = transpose(repelem(1:N_div,N_div));

% delete additional particles created by ceil() function
x_0 = x_0(1:end-(length(x_0)-p),:);
y_0 = y_0(1:end-(length(y_0)-p),:);
% move all particles slightly, to ensure they are within the box
x_0 = L.*(1/(max(x_0)) .* x_0) - 4*radius;
y_0 = L.*(1/(max(y_0)) .* y_0) - 4*radius;

%initial particle velocities
sd = sqrt((k*T_i)/(m));  % standard deviation of particles
vx_j = randn(p,sims)*sd;   % initial x vecocity in m/s
vy_j = randn(p,sims)*sd;   % initial y vecocity in m/s

Tau = min(abs(min(L./min(vx_j,vy_j))));   % time for particle to traverse box
tmax = Tau*5;  % total simulation time
tres = sigma*sqrt(m/epsilon)*0.005;  % time resolution
N_time = 2000;  % number of time steps in simulation

vx = vx_j;  % initialise velocities
vy = vy_j;  % initialise velocities

% main simulation loop
% positions, velocities, and accelerations etc. are updated according to
% the 'Verlet method'
for circ = 1:sims
    % INITIALISE VARIABLES
    
    % x- and y- positions of the 'p' particles
    x = zeros(p,N_time); 
    y = zeros(p,N_time);
    
    % set x- and y- position of all particles to their initial positions x_0
    % and y_0
    for i = 1:p
        x(i,1) = x_0(i,circ);
        y(i,1) = y_0(i,circ);
    end

    % particle accelerations
    ax = zeros(p,N_time);
    ay = zeros(p,N_time);
    % the velocity for v(t + 0.5t), as required for the Verlet Method
    vx_half = zeros(p,N_time);
    vy_half = zeros(p,N_time);
    % velocity v(t) and v(t+dt)
    vxx = zeros(p,N_time);
    vyy = zeros(p,N_time);
    % force exerted on each particle/s
    force_x = zeros(p,N_time);
    force_y = zeros(p,N_time);
    
    vxx(:,1) = vx_j(:,1);
    vyy(:,1) = vy_j(:,1);
    
    % calculation of the position of particles for t>0
    for i = 1:(N_time-1) 
        x(:,i+1) = x(:,i) + vxx(:,i).*(tres) + (1/2).*ax(:,i).* (tres).^2;
        y(:,i+1) = y(:,i) + vyy(:,i).*(tres) + (1/2).*ay(:,i).* (tres).^2;
        vx_half(:,i+1) = vxx(:,i) + (1/2).*ax(:,i).*(tres);
        vy_half(:,i+1) = vyy(:,i) + (1/2).*ay(:,i).*(tres);
          
        % distance between particle pairs
        X = x(:,i+1) - transpose(x(:,i+1));
        Y = y(:,i+1) - transpose(y(:,i+1));
        % resets for every loop i
        force_on_i_x = zeros(p,p);
        force_on_i_y = zeros(p,p);
        
        % calculate forces
        for a = 1:p
            for b = 1:p
               if a == b  % a particle exerts no force on itself
                   force_on_i_x(a,b) = 0;
                   force_on_i_y(a,b) = 0;
               else
                   x_ab = X(a,b);
                   y_ab = Y(a,b);
    
                   % pre-defined distance
                   % if two particles are within this distance, 
                   % their interaction is non-negligable, and their interaction forces must be calculated
                   r2 = x_ab.^2 + y_ab.^2;  
                   if(r2 < 3*(A/12)^(1/12))  % if within this distance, calculate forces
                       f_ab_x = (12 .* A .* x_ab)./((r2).^7) - (6 .* B .* x_ab)./((r2).^4);
                       f_ab_y = (12 .* A .* y_ab)./((r2).^7) - (6 .* B .* y_ab)./((r2).^4);
                   else  % if outside distance, force can be approximated as zero
                       f_ab_x = 0;
                       f_ab_y = 0;
                   end

                   % find force on ab'th particle
                   force_on_i_x(a,b) = force_on_i_x(a,b) + f_ab_x;  
                   force_on_i_y(a,b) = force_on_i_y(a,b) + f_ab_y; 
               end
            end
        end
        
        % sum up forces, and compute accelerations
        for particle = 1:p
            force_x(particle,i) = sum(force_on_i_x(:,particle));
            force_y(particle,i) = sum(force_on_i_y(:,particle));
            ax(particle,i+1) = -force_x(particle,i)./m;
            ay(particle,i+1) = -force_y(particle,i)./m;
        end
    
        % update vxx and vyy
        vxx(:,i+1) = vx_half(:,i+1) + (1/2) .* ax(:,i+1) .* (tres);
        vyy(:,i+1) = vy_half(:,i+1) + (1/2) .* ay(:,i+1) .* (tres);
        
        %  check if outside box, and reverse direction if so
        ind = find(x(:,i+1)>L);  
        vxx(ind,i+1) = -vxx(ind,i+1);       
        ind = find(x(:,i+1)<0);
        vxx(ind,i+1) = -vxx(ind,i+1);
        ind = find(y(:,i+1)>L);
        vyy(ind,i+1) = -vyy(ind,i+1);
        ind = find(y(:,i+1)<0);
        vyy(ind,i+1) = -vyy(ind,i+1);
        
        %  'cool' particles by reducing their speed by a factor of 0.9999
        vxx(:,i+1) = 0.9999 .* vxx(:,i+1); 
        vyy(:,i+1) = 0.9999 .* vyy(:,i+1); 
    end 
end 

%PLOTTING
hold on
axis equal
t = linspace(0,2*pi);  % full circle revolution, used for 'fill' below
% plot the initial position of the 'p' particles
for h = 1:p
    fill(x_0(h,circ)+0.5*radius*cos(t), y_0(h,circ)+0.5*radius*sin(t), [1 0 0]);
end
axis([0 L 0 L])
f = getframe;

[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,20) = 0;  % save the position as the first frame of a movie

% plot the position of the 'p' particles for t>0
for i = 1:N_time
    if round(i/250) == i/250  % update particle positions every 250th cycle of the loop, so the code isn't slowed down by constant plotting
        cla
        hold on
        axis equal
        for h = 1:p
            fill(x(h,i)+0.5*radius*cos(t),y(h,i)+0.5*radius*sin(t),[1 0 0]);
        end
        axis([0 L 0 L])
        f = getframe;
        im(:,:,1,i/250) = rgb2ind(f.cdata,map,'nodither'); %save the time steps of the plot as a movie
    end
end

imwrite(im,map,'CoolingDownParticles.gif','DelayTime',0,'LoopCount',inf) %save the movie as a gif
