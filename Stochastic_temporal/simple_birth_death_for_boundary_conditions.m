divRate = 0.25;
probSym = 0.22;
prolifCoeff = 0.0275;  % divRate * 1/2 * probSym = probProlif
diffCoeff = 0.0006875; % parameter fit by using equilibrium cell density of 40 per mm^2 (0.0275/40)
asymRate=(divRate*(1-probSym)); % rate of asymmetric division
X = 40; % stem cell population for this model equilibrium population is equal to prolif rate / diff rate
runs = 1000;
N=5000; %number of jumps
extinctArray = zeros(runs,1);  % array to contain extinction time
popsAtTimepoint = zeros(N, runs);
timepointOfFinalN = zeros(runs,1);

for k=1:runs
    t=zeros(N,1);  % array containing time
    dt=zeros(N,1); % array containing time intervals
    z=zeros(N,1);  % state variable array (contains population)
    z(1,:)= X;    % initial values input into first array of state variabl array
    for i=2:N
        if z(i-1,1) == 0
            fprintf("Extincition occurred at %d\n", t(i-1));
            extinctArray(k) = t(i-1);
            break
        else
        x=z(i-1,1); % set pop of stem cells as previous size size
        rates=[prolifCoeff*x asymRate*x diffCoeff*x*x]; %set birth and death rates based on current pop
        R=[1 0 -1]; % set change in pop size for each of the two rates
        lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the total rate of all possible events
        dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
        t(i)=t(i-1)+dt(i-1); %update time array
        rates=rates/lam; %normalise the rates
        reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
        z(i,:)=z(i-1,:)+R(:,reac)'; %update pop size based on event that occured
        %pop of X at step i = pop at step i-1 + stoichiometry matrix
        popsAtTimepoint(round(t(i))+1, k) = x; %%this might be an issue, if t is larger than N indexing will be out of bounds
        end
        if i == N
            timepointOfFinalN(k) = t(i);
        end
    end
    plot(t, z);
    title('Population of Stem Cells in 1mm^2 patches of skin')
    xlabel('Time (weeks)')
    ylabel('Stem cell population')
    hold on
end

disp(popsAtTimepoint(1:2, :));
disp(timepointOfFinalN);
finalTime = round(min(timepointOfFinalN));
trimmedPopsAtTimepoint = popsAtTimepoint(1:finalTime,:);
csvwrite("populations_for_PDE_boundary"  , trimmedPopsAtTimepoint);
disp(width(trimmedPopsAtTimepoint))
