clear all;
close all;
clc;

% print out test identifer so can track while running multiple tests
trialparams.test = 1
trialparams.trialnum = 0; %1;

trialparams.outputdir = strcat("test_", num2str(trialparams.test));
if exist(trialparams.outputdir, 'dir')
    rmdir(trialparams.outputdir,'s');
end
mkdir(trialparams.outputdir);

% reproduce sine oscilators from Supplementary Figure 2
if trialparams.test == 0
    trialparams.test_name = "oscilators";
    trialparams.test_title = "sine oscilators";
    
    % values from Supplementary Table 1
    trialparams.omegaparams = [
        [5*10^3, 5*10^3]; % sine
        [5*10^3, 4*10^3]; % sawtooth
        [1*10^4, 9*10^3]; % sine product
        [1*10^4, 8*10^3]; % noisy sine product
        [1*10^4, 1*10^4]; % oscillator 1
        [1*10^4, 1*10^4] % oscillator 2
    ];

    for i=1:6
        trialparams.supervisor_choice = i;
        trialparams.G = trialparams.omegaparams(i, 1);
        trialparams.Q = trialparams.omegaparams(i, 2);
        
        svparam.normalize = 0;
        
        if i == 5 || i == 6
            % normalize the oscilators since amplitude is large causes 
            svparam.normalize = 1;
            
            % train the oscilators longer
            trialparams.trainlength = 25;
        end

        IZFORCESINE
        
        trialparams.trialnum = trialparams.trialnum + 1;
    end
% Check training robustness across trails looking for random effects
elseif trialparams.test == 1
    trialparams.test_name = "robustness";
    trialparams.test_title = "sine";
    
    trialparams.numtrials = 30; %10;

    for i=1:trialparams.numtrials
        IZFORCESINE
        trialparams.trialnum = trialparams.trialnum + 1;
    end
% elseif trialparams.test == -1 % temporary test noticed that the saw tooth wasn't reproducing well
%     trialparams.supervisor_choice=2;
%     trialparams.Q=4*10^3;
%     trialparams.test_name = "robustness";
%     trialparams.test_title = "sawtooth";
%     trialparams.numtrials = 10;
% 
%     for i=1:trialparams.numtrials
%         IZFORCESINE
%         trialparams.trialnum = trialparams.trialnum + 1;
%     end
% elseif trialparams.test == -2 % curious about long term behavoir
%     trialparams.supervisor_choice=2;
%     trialparams.Q=4*10^3;
%     trialparams.test_name = "robustness";
%     trialparams.test_title = "sawtooth long";
%     trialparams.numtrials = 10;
%     trialparams.simlength = 15;
% 
%     for i=1:trialparams.numtrials
%         IZFORCESINE
%         trialparams.trialnum = trialparams.trialnum + 1;
%     end

% vary signal
elseif trialparams.test > 1
    trialparams.test_name = "supervisors";

    if trialparams.test == 2
        trialparams.test_title = "supervisors amplitude";
        trialparams.rng = logspace(-2,2,5); %[10^-2:10^2]
    elseif trialparams.test == 3
        trialparams.test_title = "supervisors frequency";
        trialparams.rng = 5*logspace(-2,2,5);
    elseif trialparams.test == 4
        trialparams.test_title = "supervisors phase shift";
        trialparams.rng = linspace(0,1,5);
    elseif trialparams.test == 5
        trialparams.test_title = "supervisors vert shift";
        trialparams.rng = logspace(-2,2,5);
    elseif trialparams.test == 6
        trialparams.test_title = "supervisors complexity";
        trialparams.rng = [1:1:6];
    end

    for i=1:length(trialparams.rng)
        if trialparams.test == 2
            svparam.amplitude = trialparams.rng(i);
        elseif trialparams.test == 3
            svparam.frequency = trialparams.rng(i);
        elseif trialparams.test == 4
            svparam.phaseshift = trialparams.rng(i);
        elseif trialparams.test == 5
            svparam.verticalshift = trialparams.rng(i);
        elseif trialparams.test == 6
            trialparams.supervisor_choice = trialparams.rng(i);
            svparam.normalize = 0;
        
            % normalize the oscilators since amplitude is large causes 
            if i == 5 || i == 6
                svparam.normalize = 1;
            end
        end

        IZFORCESINE
        trialparams.trialnum = trialparams.trialnum + 1;
    end
end

% If there are more than 9 trials, the files are not sorted in in numerical
% order so sort so that trials are aligned.
trials=dir(strcat(trialparams.outputdir,'/*_izforcesine.mat'));

names = {trials.name};
str  = sprintf('%s#', names{:});
num = sscanf(str, '%d_izforcesine.mat#');
[~, index] = sort(num);
trials = names(index);

for i=1:length(trials)
    load(strcat(trialparams.outputdir, "/",trials{i}),'trialinfo','simulationinfo');
    NRMSE(i)=trialinfo.NRMSE;
    ET(i)=trialinfo.ElapsedTime;
    AFR(i)=trialinfo.AverageFiringRate;
end

hyperparameters = strcat("{\fontsize{8} N: ", num2str(simulationinfo.N), ", G: ", num2str(simulationinfo.G), ", Q: ", num2str(simulationinfo.Q), "}");
duration = strcat("{\fontsize{8}Simulation Duration: ", num2str(simulationinfo.length), " s, Sample time: ", num2str(simulationinfo.sampletime), " ms, Time steps: ", num2str(simulationinfo.timesteps), "}");
training = strcat("{\it\fontsize{8}[Pre-training: ", num2str(simulationinfo.trainstart), " s, Training: ", num2str(simulationinfo.trainstop-simulationinfo.trainstart), " s, Post-Training: ", num2str(simulationinfo.length-simulationinfo.trainstop), " s]}");

fig=figure;
bar(NRMSE);
yline(mean(NRMSE),'--',{'mean';num2str(mean(NRMSE))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " NRMSE across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('NRMSE');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_rmse"),'png');

fig=figure;
bar(NRMSE);
set(gca, 'YScale', 'log');
yline(mean(NRMSE),'--',{'mean';num2str(mean(NRMSE))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " log(NRMSE) across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('NRMSE');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_logrmse"),'png');

fig=figure;
bar(ET);
yline(mean(ET),'--',{'mean';num2str(mean(ET))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " Elapsed Time across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('Time (s)');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_et"),'png');

fig=figure;
bar(ET);
set(gca, 'YScale', 'log');
yline(mean(ET),'--',{'mean';num2str(mean(ET))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " log(Elapsed Time) across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('Time (s)');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_loget"),'png');

fig=figure;
bar(AFR);
yline(mean(AFR),'--',{'mean';num2str(mean(AFR))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " Average Firing Rate across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('rate');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_afr"),'png');

fig=figure;
bar(AFR);
set(gca, 'YScale', 'log');
yline(mean(AFR),'--',{'mean';num2str(mean(AFR))},'LabelVerticalAlignment','middle');
title({strcat("{\bf", trialparams.test_title, " log(Average Firing Rate) across trials}");hyperparameters;duration;training},'FontWeight','Normal');
xlabel('trial');
ylabel('rate');
saveas(fig,strcat(trialparams.outputdir, "/", trialparams.test_name, "_logafr"),'png');
