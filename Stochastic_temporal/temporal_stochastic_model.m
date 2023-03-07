% Initialise run variables
X=100;  % initial conditions
k=60; % birth rate coefficient
lambda=0.1; % death rate coefficient

N=1000; % Number of jumps
t=zeros(N,1);  % array containing time
dt=zeros(N,1); % array containing time intervals
z=zeros(N,1);  % state variable array (contains population)
z(1,:)=[X];    % initial values input into first array of state variabl array

for i=2:N
    x=z(i-1,1); % x in the previous state variable
    rates=[k lambda*x]; %the parameters used to work out rates
    R=[1 -1]; % no idea what this is
    lam=sum(rates); % 
    dt(i-1)= -log(rand)/lam;
    t(i)=t(i-1)+dt(i-1);
    rates=rates/lam;
    reac=1+sum(rand>cumsum(rates));
    z(i,:)=z(i-1,:)+R(:,reac)';
end

plot(t,z);

% Generate weighted mean (last element not used)
transpose(dt)*z/sum(dt)

% Generate weighted sum of squares (last element not used)
var(z,dt)


%% implement a birth death process with a second reaction where it goes back to itself
%% plot a graph of how this p(asymmetric division) affects the variance
%% set birthrate to 1 to fit death rate parameter
