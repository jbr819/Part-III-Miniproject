% Initialise run variables
prolifRate=0.0275; % birth rate coefficient
diffRate=0.0275; % death rate coefficient
asymRate=0.195; % rate of asymmetric division
keratinLoss=0;
X = [10, 0]; % for this model equilibrium population is equal to prolif rate / diff rate

N=100000; % Number of jumps
t=zeros(N,1);  % array containing time
dt=zeros(N,1); % array containing time intervals
z=zeros(N,2);  % state variable array (contains population)
z(1,:)=[X];    % initial values input into first array of state variabl array

for i=2:N
    x=z(i-1,1); % set previous pop of stem cells as current size
    y = z(i-1, 2); % set previous pop of keratinocyte cells as current size
    rates=[prolifRate*x asymRate*x diffRate*x]; %set birth and death rates based on current pop
    R=[1,-1; 0, 1; -1, 2]; % set change in pop size for each of the two rates
    lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the totla rate of all possible events
    dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
    t(i)=t(i-1)+dt(i-1); %update time array
    rates=rates/lam; %normalise the rates
    reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
    z(i,:)=z(i-1,:)+R(:,reac)'; %update pop size based on event that occured
end

plot(t,z);
hold on

% Generate weighted mean (last element not used)
transpose(dt)*z/sum(dt)

% Generate weighted sum of squares (last element not used)
var(z,dt)


%% implement a birth death process with a second reaction where it goes back to itself
%% plot a graph of how this p(asymmetric division) affects the variance and the extinction distribution
%% set birthrate to 1 to fit death rate parameter
