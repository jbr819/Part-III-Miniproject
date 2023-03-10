% Initialise run variables
divRate = 0.25;
probSym = 0.22; 
prolifDiffRatio = 0.5;
prolifRate=(divRate*probSym*(1-prolifDiffRatio)); % birth rate coefficient
diffRate=(divRate*probSym*prolifDiffRatio); % death rate coefficient
asymRate=(divRate*(1-probSym)); % rate of asymmetric division
X = 400000; % for this model equilibrium population is equal to prolif rate / diff rate
runs = 1;
extinctArray = zeros(runs,1);  % array to contain extinction time

for k=1:runs
    N=1000000; % Number of jumps
    t=zeros(N,1);  % array containing time
    dt=zeros(N,1); % array containing time intervals
    z=zeros(N,2);  % state variable array (contains population)
    z(2,:)=X;    % initial values input into first array of state variabl array
    for i=2:N
        if z(i-1) == 0
            fprintf("Extincition occurred at %d\n", t(i-1));
            extinctArray(k) = t(i-1);
            break
        else
        x=z(i-1,1); % set previous pop of stem cells as current size
        rates=[prolifRate*x asymRate*x diffRate*x*x]; %set birth and death rates based on current pop
        R=[1 0 -1]; % set change in pop size for each of the two rates
        lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the total rate of all possible events
        dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
        t(i)=t(i-1)+dt(i-1); %update time array
        rates=rates/lam; %normalise the rates
        reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
        z(i,:)=z(i-1,:)+R(:,reac)'; %update pop size based on event that occured
        end
    end
    plot(t, z);
    hold on
end

disp(extinctArray);
histogram(extinctArray, 'Normalization','probability')



% Generate weighted mean (last element not used)
% transpose(dt)*z/sum(dt)

% Generate weighted sum of squares (last element not used)
% var(z,dt)
