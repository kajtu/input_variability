function [simAmi,sol] = simulate_avatarbelen(theta,constants,options,numHeartBeats,ind,simtime)
% Do a simulation of the avatar model. First, set parameter values, then
% simulate a number of heartbeats one at a time.
%step = 0.001; %good resolution when plotting simulation
%step = T/39; %comparing to data (since there are 40 timeframes)
%ind is a struct containing indexes for the parameters
%options contains simulation options
T = constants(ind.T);

% First simulation of one heartbeat - to get the sizes of x and y
simtime = simtime+T;
options.tstart = T;
sol = simulate_avatar_belen(simtime,theta, constants, [], options);

options.x0 = sol.x(end,:)';
simlen = length(sol.t)-1;
simAmi.x = zeros(simlen*numHeartBeats,size(sol.x,2));
simAmi.y = zeros(simlen*numHeartBeats,size(sol.y,2));
simAmi.t = zeros(simlen*numHeartBeats,size(sol.t,2));

simAmi.x(1:length(sol.t),:) = sol.x;
simAmi.y(1:length(sol.t),:) = sol.y;
simAmi.t(1:length(sol.t)) = sol.t;

% The rest of the simulations
for i = 1:numHeartBeats-1
    sol = simulate_avatar_belen(simtime,theta, constants, [], options);

    %start next sim with end values from this sim
    options.x0 = sol.x(end,:)';
    
    %save this sim together with the rest
    simAmi.x(simlen*i+2:simlen*(i+1)+1,:) = sol.x(2:end,:);
    simAmi.y(simlen*i+2:simlen*(i+1)+1,:) = sol.y(2:end,:);
    simAmi.t(simlen*i+2:simlen*(i+1)+1)   = sol.t(2:end)+T*i;
end

 % To compare with data we use the last heartbeat when a "steady state" is established.
 % We need to round the time to 3 digits, and start at t=0 to be able to
 % compare with data.
sol.t = round(sol.t-T,3);

end