divRate = 0.25;
probSym = 0.22;
prolifCoeff = 0.0275;
diffCoeff = 0.0006875;
flakingCoeff = 0.25; % each keratinocyte has a 
asymRate=(divRate*(1-probSym)); % rate of asymmetric division
X = [40, 200]; % stem cell population for this model equilibrium population is equal to prolif rate / diff rate
K = 200; % the keratinocytpe population
runs = 5;
extinctArray = zeros(runs,1);  % array to contain extinction time

for k=1:runs
    N=1000; % Number of jumps
    t=zeros(N,1);  % array containing time
    dt=zeros(N,1); % array containing time intervals
    z=zeros(N,2);  % state variable array (contains population)
    z(2,:)= X;    % initial values input into first array of state variabl array
    for i=2:N
        if z(i-1,1) == 0
            fprintf("Extincition occurred at %d\n", t(i-1));
            extinctArray(k) = t(i-1);
            break
        else
        x=z(i-1,1); % set pop of stem cells as previous size size
        y=z(i-1,2); % set pop of keratinocytes to previous pop
        rates=[prolifRate*x asymRate*x diffRate*x*x flakingCoeff*y]; %set birth and death rates based on current pop
        R=[1,0; 0,1; -1,2; 0,-1]; % set change in pop size for each of the two rates
        lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the total rate of all possible events
        dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
        t(i)=t(i-1)+dt(i-1); %update time array
        rates=rates/lam; %normalise the rates
        reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
        z(i,1)=z(i-1,1)+R(:,reac)'; %update pop size based on event that occured
        z(i,2)=z(i-1,2)+R(:,reac)'; %update pop size based on event that occured
        end
    end
    plot(t, z);
    hold on
end

%disp(extinctArray);
%histogram(extinctArray, 'Normalization','probability')
