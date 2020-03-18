%% Force Method with Izhikevich Network 
clearvars -except trialparams svparam;
close all
%clc 

% keep trials separate
if exist('trialparams','var') && isfield(trialparams,'trialnum')
    trialnum = trialparams.trialnum;
else
    trialnum = 0;
end

trialinfo.trialnum = trialnum;

% what to plot
plotvars.timeseries = 1;
plotvars.decoders = 1;
plotvars.raster = 1;
plotvars.voltage = 1;
plotvars.eigenvalues = 1;
plotvars.rmse = 1;
plotvars.xcorr = 1;

% what to calculate
statvars.afr = 1;
statvars.nrmse = 1;
statvars.xcorr = 1;
statvars.meanfreq = 1;

%% simulation parameters
%rng(1); %set random number generator seed
trialinfo.RandomNumberCheck = rand();

% time in seconds
if exist('trialparams','var') && isfield(trialparams,'trainstart')
    simulationinfo.trainstart = trialparams.trainstart;
else
    % default pre training length 5 s
    simulationinfo.trainstart = 5;
end

if exist('trialparams','var') && isfield(trialparams,'trainlength')
    simulationinfo.trainstop = simulationinfo.trainstart + trialparams.trainlength;
else
    % default training length 5 s
    simulationinfo.trainstop = simulationinfo.trainstart + 5; 
end

if exist('trialparams','var') && isfield(trialparams,'simlength')
    simulationinfo.length = simulationinfo.trainstop + trialparams.simlength;
else
     % default simulation length is 5 s
    simulationinfo.length = simulationinfo.trainstop + 5;
end

T = simulationinfo.length*1000; % simulate for 15 seconds (total time T in ms)
dt = 0.04; % convert ms to times steps (e.g. 4 ms / 100 time steps)
nt = round(T/dt); % calculate number of time steps for entier simulation

simulationinfo.sampletime = dt; % in ms
simulationinfo.timesteps = nt;

simulationinfo.N = 2000;  % number of neurons 

%% training paramters
% start traiing at 5 seconds and end at 10 seconds

imin = round(simulationinfo.trainstart*1000/dt); %time before starting RLS, gets the network to chaotic attractor 
icrit = round(simulationinfo.trainstop*1000/dt); %end simulation at this time step 

%% izhikevich parameters
load IZPARAMS

%% voltage initailization
v = vr+(vpeak-vr)*rand(simulationinfo.N,1);
v_ = v; % store previous voltage

% adaptation current for izhikevich neurons
u = zeros(simulationinfo.N,1);

%% variables for synapse integration  

% static term
Ipsc = zeros(simulationinfo.N,1);
hpsc = zeros(simulationinfo.N,1); 
staticTerm = zeros(simulationinfo.N,1);

% feedback term
r = zeros(simulationinfo.N,1);
hr = zeros(simulationinfo.N,1);

%% supervisor signal
% allow calling script to dictate supervisor choice
if exist('trialparams','var') && isfield(trialparams,'supervisor_choice')
    svparam.choice = trialparams.supervisor_choice;
else
    svparam.choice = 1; % default
end
% svparam.amplitude = 1;
% svparam.freqence = 5;
% svparam.phaseshift = 0;
% svparam.verticalshift = 0;
svparam.duration = nt;
svparam.sampletime = dt;

supervisornames = ["sine","sawtooth","sineproduct","noisysineproduct","oscillator1","oscillator2"];
simulationinfo.supervisor = supervisornames(svparam.choice);

x = supervisor(svparam); % [m, nt]
%plot(x(1:1000/dt)); %plot 1 sec of supervisor
m = size(x,1); % dimensionality of the supervisor

%% weight parameters
%load OMEGAPARAMS
p = 0.1; %sparsity 

if exist('trialparams','var') && isfield(trialparams,'G')
    G = trialparams.G;
else
    G = 5*10^3; % default
end

if exist('trialparams','var') && isfield(trialparams,'Q')
    Q = trialparams.Q;
else
    Q = 5*10^3; % default
end

% store into simulation information
simulationinfo.G = G;
simulationinfo.Q = Q;

% sparse static weight matrix for chaotic dynamics
omega0 = (randn(simulationinfo.N,simulationinfo.N)).*(rand(simulationinfo.N,simulationinfo.N)<p)/(p*sqrt(simulationinfo.N)); % ~N(0, (p*sqrt(N))^-1) 

% feedback encoder decoders
eta = (2*rand(simulationinfo.N,m)-1);  % encoders ~U(-1,1)
phi = zeros(simulationinfo.N,m); % decoders

xhat = zeros(m,1);  %initial approximant

%% rls parameters
delt = 2/dt; % run RLS every 2ms
lambda=1/2;
P = eye(simulationinfo.N)/lambda; %initial correlation matrix, coefficient is the regularization constant as well 

%% variables for plotting
output = zeros(nt,m);  %store the approximant 

if statvars.nrmse
    SE = zeros(nt,m); % squared error
end

if plotvars.rmse
    RMSE = zeros(nt,m); % root mean squared error
end

if plotvars.decoders
    numsamples = 5;
    decoders = zeros(nt,numsamples); %store the decoders 
end

if plotvars.voltage 
    numsamples = 5;
    voltagetrace = zeros(nt,numsamples); %store voltage
    adaptationtrace = zeros(nt,numsamples); %store adaptation
end

if plotvars.raster || statvars.afr
    % size 5*nt is arbitarily large
    tspike = zeros(5*nt,2);  %If you want to store spike times, 
    totnumspikes = 0; %count total number of spikes
end

%% SIMULATION
tic;

for i = 1:1:nt 
    %% EULER INTEGRATE
    I = Ibias + Ipsc + Q*eta*xhat;  %postsynaptic current 
    v = v + dt*(( k.*(v-vr).*(v-vt) - u + I))/C ; % v(t) = v(t-1)+dt*v'(t-1)
    u = u + dt*(a*(b*(v_-vr)-u)); %same with u, the v_ term makes it so that the integration of u uses v(t-1), instead of the updated v(t)

    %% double exponential synapse 
    
    % split synapses into calculations for static and feedback terms
    
    % static term contributions from spiking neurons
    index = find(v>=vpeak);

    if ~isempty(index)
        % compute static increase in current from spiking neurons
        staticTerm = G*sum(omega0(:,index),2);
        
        if plotvars.raster || statvars.afr
            % track all spikes at this time step by index and time (ms) of spike
            tspike(totnumspikes+1:totnumspikes+length(index),:) = [index,0*index+dt*i];
            totnumspikes = totnumspikes + length(index); 
        end
    end
    
    Ipsc = Ipsc*exp(-dt/tr) + hpsc*dt;
    hpsc = hpsc*exp(-dt/td) + staticTerm*(~isempty(index))/(tr*td);

    % feedback term contributions from spiking neurons
    r = r*exp(-dt/tr) + hr*dt; 
    hr = hr*exp(-dt/td) + (v>=vpeak)/(tr*td);
    xhat = phi'*r; % output 
    
    % error signal for output used to learn decoders phi
    e = xhat - x(:,i); % error 
    
    %% RLS 
    if mod(i,delt)==1 && i > imin && i < icrit 
        g = P*r; % gain vector
        phi = phi - (g*e'); % update the decoder
        P = P - ((g)*(g'))/( 1 + (r')*(g)); % update the inverse correlation matrix
    end

    %% reset
    u = u + d*(v>=vpeak);  %implements set u to u+d if v>vpeak, component by component. 
    v = v+(vreset-v).*(v>=vpeak); %implements v = c if v>vpeak add 0 if false, add c-v if true, v+c-v = c
    v_ = v;  % sets v(t-1) = v for the next itteration of loop

    %% Store 
    output(i,:) = xhat';
    
    if statvars.nrmse
        SE(i,:) = e.^2;
    end
    
    if plotvars.decoders
        decoders(i,:) = phi(1:numsamples);
    end
    
    if plotvars.voltage
        voltagetrace(i,:) = v(1:numsamples)';
        adaptationtrace(i,:) = u(1:numsamples)';
    end
    
    if plotvars.rmse
        RMSE(i,:) = rms(e); %sqrt(mean(e.^2));
    end
end

trialinfo.ElapsedTime=toc;

%% Average firing rate of neurons

if statvars.afr
    % remove empty rows from the tspike matrix
    tspike = tspike(tspike(:,2)~=0,:); 
    % find the post training spikes
    M = tspike(tspike(:,2)>dt*icrit); 
    % average of spikes per neuron per second post training
    trialinfo.AverageFiringRate = length(M)/simulationinfo.N*1000/(T-dt*icrit);
end

if statvars.nrmse
    % select only the post training part of supervisor
    replayx=x(icrit:nt);
    
    meanx = mean(replayx);
    
    % retify wave to use only the half cycle to compute average
    if meanx == 0
        meanx = mean(replayx(replayx>0));
    end
        
    % compute NRMSE
    trialinfo.NRMSE=1/meanx*sqrt(1/(nt-icrit)*trapz(SE(icrit:nt)));
end

if statvars.xcorr
    trialinfo.corr=xcorr(output(icrit:end)',x(icrit:end),0,'coeff');
end

if statvars.meanfreq
    trialinfo.meanfreq=meanfreq(output(icrit:end)',25000);
end

%% Plotting

if exist('trialparams','var') && isfield(trialparams,'outputdir')
    trialdir = trialparams.outputdir;
else
    trialdir = "trials";
    
    if ~exist(trialdir, 'dir')
        mkdir(trialdir);
    end
end

if plotvars.timeseries
    fig=figure('Name','1','Visible', 'off');
    %plot((1:1:nt)*dt/1000,x(:,1:1:nt),'k','LineWidth',2), hold on
    plot((icrit:1:nt)*dt/1000,x(:,icrit:1:nt),'k','LineWidth',2), hold on
    %plot((1:1:nt)*dt/1000,output(1:1:nt,:),'b--','LineWidth',2), hold off
    plot((icrit:1:nt)*dt/1000,output(icrit:1:nt,:),'b--','LineWidth',2), hold off
    xlabel('Time (s)')
    ylabel('Amplitude')
    legend('x','$\hat{x}(t)$','Interpreter','LaTeX')
    %xlim([icrit*dt,nt*dt]/1000) % testing
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_timeseries'),'png');
end

if plotvars.decoders
    fig=figure('Name','2','Visible', 'off');
    plot((1:1:nt)*dt/1000,decoders(1:1:nt,:));
    %xline(imin,'--',{'training start'});
    %xline(icrit,'--',{'training stop'});
    set(gca,'YTick', []);
    ylabel('Decoders');
    xlabel('Time (s)');
    title(['First ', num2str(numsamples) ,' Decoders \phi']);
    % xlim([imin*dt,icrit*dt]/1000) % training
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_decoders'),'png');
end

if plotvars.raster
    rasterLim = 100;
    fig = figure('Name','3','Visible', 'off');
    plot(tspike(1:totnumspikes,2),tspike(1:totnumspikes,1),'k.');
    set(gca, 'XTick', get(gca, 'XTick'), 'XTickLabel', get(gca, 'XTick')/1000);
    ylim([0,rasterLim]);
    xlabel('Time (s)');
    ylabel('Neuron index');
    title(['Raster Plot Neurons 1-', num2str(rasterLim)]);
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_raster'),'png');
end

if plotvars.voltage
    fig = figure('Name','4','Visible', 'off');
    for j = 1:1:numsamples
        plot((1:1:nt)*dt/1000,voltagetrace(1:1:nt,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([icrit*dt,nt*dt]/1000); % testing
    xlabel('Time (s)');
    ylabel('Neuron Index'); 
    title('post learning voltage trace');
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_postvoltage'),'png');

    fig = figure('Name','5','Visible', 'off');
    for j = 1:1:numsamples
        plot((1:1:nt)*dt/1000,voltagetrace(1:1:nt,j)/(vpeak-vreset)+j), hold on 
    end
    xlim([0,imin*dt]/1000); % pre training
    xlabel('Time (s)');
    ylabel('Neuron Index');
    title('pre-learning voltage trace');
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_prevoltage'),'png');
end

if plotvars.eigenvalues
    fig = figure('Name','6','Visible', 'off');
    Z = eig(G*omega0+Q*eta*phi'); %eigenvalues after learning 
    Z2 = eig(G*omega0); %eigenvalues before learning 
    plot(Z2,'r.'), hold on 
    plot(Z,'k.') 
    legend('Pre-Learning','Post-Learning')
    xlabel('Re \lambda')
    ylabel('Im \lambda')
    title('eigenvalue spectrum')
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_eigenspectrum'),'png');
end

if plotvars.rmse
    fig = figure('Name','6','Visible', 'off');
    plot((1:1:nt)*dt/1000,RMSE)
    title('RMSE')
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_rmse'),'png');
end

if plotvars.xcorr
    fig = figure('Name','7','Visible', 'off');
    plot(x(icrit:end),output(icrit:end)');
    title('correlation');
    xlabel('supervisor');
    ylabel('generated');
    saveas(fig,strcat(trialdir, '/',num2str(trialnum), '_corr'),'png');
end

% print out information
simulationinfo
trialinfo

% SAVE variables
save(strcat(trialdir, '/', num2str(trialnum),'_izforcesine.mat'), 'simulationinfo', 'trialinfo');
save(strcat(trialdir, '/', num2str(trialnum),'_output.mat'), 'output');