divRate = 0.25; %(cell week-1)
probSym = 0.22; % lambda
ProbProlif = 0.7;
ProbDiff = 1-ProbProlif;
eqmStemDensity = 40;
prolifRateCoeff = divRate * probSym * ProbProlif; % divRate * 1/2 * probSym = probProlif
diffRateCoeff = divRate * probSym * ProbDiff;  % this is too well coded, needs to allow eqm stem cell density to vary with different proliferation rate constants
asymRateCoeff =(divRate*(1-probSym)); % rate of asymmetric division
competitiveDiff = 0.00055; % a rate constant controlling rate of winner-loser to achieve ss pop of 40 (per mm^2)
keratinShedRateCoeff = 0.065; % (cells per week per mm2)
X = [40, 2700]; % stem cell = index 1, keratinocyte = index 2
runs = 1;
N=50000; %number of jumps
extinctArray = zeros(runs,1);  % array to contain extinction time

startTimeOfLesion = 50;
endTimeOfLesion = 70;

for k=1:runs
    t=zeros(N,1);  % array containing time
    dt=zeros(N,1); % array containing time intervals
    z=zeros(N,2);  % state variable array (contains population)
    z(1,:)= X;     % initial values input into first array of state variable array
    for i=2:N
        if z(i-1,1) == 0
            fprintf("Extincition occurred at %d\n", t(i-1));
            extinctArray(k) = t(i-1);
            break
        else
            if (t(i-1) > startTimeOfLesion) && (t(i-1) < endTimeOfLesion)
                disp("Lesioning\n");
                ProbProlif = 0.3;
                ProbDiff = 1-ProbProlif;
                eqmStemDensity = 40;
                prolifRateCoeff = divRate * probSym * ProbProlif; % divRate * 1/2 * probSym = probProlif    
                diffRateCoeff = divRate * probSym * ProbDiff;  % this is too well coded, needs to allow eqm stem cell density to vary with different proliferation rate constants
                asymRateCoeff =(divRate*(1-probSym)); % rate of asymmetric division

                x=z(i-1,:); % set pop of stem/diff cells as previous size size
                rates=[prolifRateCoeff*x(1) asymRateCoeff*x(1) diffRateCoeff*x(1) keratinShedRateCoeff*x(2) competitiveDiff*x(1)*x(1)]; %set birth and death rates based on current pop
                R=[1,0; 0,16; -1,32; 0,-1; -1,32]; % set change in pop size for each of the two rates
                lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the total rate of all possible events
                dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
                t(i)=t(i-1)+dt(i-1); %update time array
                rates=rates/lam; %normalise the rates
                reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
                z(i,2)=z(i-1,2)+R(reac,2)'; %update pop size based on event that occured
                z(i,1)=z(i-1,1)+R(reac,1)'; %update pop size based on event that occured
            else
                ProbProlif = 0.7;
                ProbDiff = 1-ProbProlif;
                eqmStemDensity = 40;
                prolifRateCoeff = divRate * probSym * ProbProlif; % divRate * 1/2 * probSym = probProlif    
                diffRateCoeff = divRate * probSym * ProbDiff;  % this is too well coded, needs to allow eqm stem cell density to vary with different proliferation rate constants
                asymRateCoeff =(divRate*(1-probSym)); % rate of asymmetric division
                

                x=z(i-1,:); % set pop of stem/diff cells as previous size size
                rates=[prolifRateCoeff*x(1) asymRateCoeff*x(1) diffRateCoeff*x(1) keratinShedRateCoeff*x(2) competitiveDiff*x(1)*x(1)]; %set birth and death rates based on current pop
                R=[1,0; 0,16; -1,32; 0,-1; -1,32]; % set change in pop size for each of the two rates
                lam=sum(rates); % calculate the total rate. Parameter of expo. dist. used to determine next event is the total rate of all possible events
                dt(i-1)= -log(rand)/lam; %set time interval based on the exponenetial distribution
                t(i)=t(i-1)+dt(i-1); %update time array
                rates=rates/lam; %normalise the rates
                reac=1+sum(rand>cumsum(rates)); %determine which event occured based on normalised rates
                z(i,2)=z(i-1,2)+R(reac,2)'; %update pop size based on event that occured
                z(i,1)=z(i-1,1)+R(reac,1)'; %update pop size based on event that occured
            end
        end
    end
    yyaxis left
    plot(t(2:N), z(2:N,1), "Color", "Blue")
    xlim([0,260])
    ylim([0,60])
    ylabel("Population of Epidermal Stem Cells", "Color","Blue");
    yyaxis right
    plot(t(2:N), z(2:N,2))
    ylim([0,3500])
    ylabel('Population of Keratinocytes')
    title('Population of cells in a healthy 1mm^2 patch of skin')
    xlabel('Time (weeks)')
end


