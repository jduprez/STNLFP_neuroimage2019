% The following code has been used to analyze the data presented in the paper:
%'Subthalamic nucleus local field potentials recordings reveal subtle effects of promised reward during conflict resolution in Parkinson's disease' (Duprez et al., XXXX)
% Time-frequency analyses are done using complex morlet wavelet convolution based on published equations (Cohen 2014 - Analyzing neural time series data).
% This code has been adapted from Mike X Cohen's one.
% Although at many places efficiency of this code can be improved, it does correct computations.


% This code applies to the data that has been epoched from -1 to 2s surrounding the onset of the imperative (simon) stimulus
% Changes in condition number will make it work for the data epoched around the onset of the reward cue.

ncond = 6; % ncond = 3 if reward cue epoching

% Soft-coded convolution parameters
min_frex  =  1; % lowest frequency
max_frex  = 40; % highest frequency
num_frex  = 50; % n frequencies

% Range for number of wavelet cycles
nCycRange = [4 10];

% Timestamps to use for visualization (downsampling), careful not to downsample before convolution !
times2save = -0.3:0.01:2;


%% Load data
% Preprocessed data for each STN is stored in a 3D matrix (bipolar montage*time*trials)

matfiles = dir(fullfile(['*.mat']));
numfiles = length(matfiles);

% Initialize output time-frequency power matrix (STN*conditions*bipolar montage*frequency*time)
% There are 6 conditions in the cognitive control epoching (3 reward * 2 congruence)
tf = zeros(numfiles, ncond, 3, num_frex, size(times2save,2));


% Initialize output time-frequency itpc matrix (STN*conditions*bipolar montage*frequency*time)
itpcz = zeros(numfiles, ncond, 3, num_frex, size(times2save,2));


for subno = 1:numfiles % loop over STNs
    
    load(matfiles(subno).name)
    
    
    % Convert to double precision (data were stored in simple precision to avoid too large data files)
    LFP.data = double(LFP.data); % Epochs
    LFP.baseline = double(LFP.baseline); % Corresponding baseline (which has been stored in the LFP structure of the cognitive control epoching but comes from the time window before the reward stimulus)
    
    % Find indices of times2save for temporal downsampling (for visualization only)
    times2saveidx = dsearchn(LFP.time',times2save');
    
    % Create logarithmically spaced frequency vector
    frex = logspace(log10(min_frex),log10(max_frex),num_frex);
    
    % wavelet parameters
    s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),num_frex)./(2*pi.*frex);
    t = -2:1/LFP.srate:2;
    halfwave = floor((length(t)-1)/2);
    
    % Data convolution parameters
    nData = LFP.trial*LFP.pnts;
    nWave = length(t);
    nConv = nData+nWave-1;
    
    % Baseline convolution parameters
    nData_base = size(LFP.baseline,2)*size(LFP.baseline,3);
    nConv_base = nData_base+nWave-1;
    
    
    % Data wavelets
    cmwX = zeros(num_frex,nConv);
    for fi=1:num_frex
        cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
        cmwX(fi,:) = cmw ./ max(cmw); % amplitude-normalize in the frequency domain
    end
    
    % Baseline wavelets
    cmwX_base = zeros(num_frex,nConv_base);
    for fi=1:num_frex
        cmw_base = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv_base);
        cmwX_base(fi,:) = cmw_base ./ max(cmw_base); % amplitude-normalize in the frequency domain
    end
    
    
    % Run convolution for each bipolar montage and for each frequency
    
    for chani = 1:size(LFP.data,1) % loop over bipolar montage
        dataX = fft( reshape(LFP.data(chani,:,:),1,nData) ,nConv); %% concatenates all the trials of the first channel
        baseX = fft( reshape(LFP.baseline (chani,:),1,nData_base) ,nConv_base);
        
        for fi=1:num_frex % loop over frequencies
            as = ifft( dataX.*cmwX(fi,:) );
            as = as(halfwave+1:end-halfwave);
            as = reshape(as,1,LFP.pnts,LFP.trial);
            
            as_base = ifft( baseX.*cmwX_base(fi,:) );
            as_base = as_base(halfwave+1:end-halfwave);
            as_base = reshape(as_base,1,size(LFP.baseline,2),size(LFP.baseline,3));
            
            % condition-average baseline
            basepow = mean(squeeze(mean(abs(as_base).^2,2)),1);
            baseitpc = mean(abs(mean(squeeze(exp(1i*angle(as_base))),2)),1);
            
            
            for condi = 1:ncond % In the LFP structure, LFP.cond is a vector containing the condition number of the trial
                tf(subno,condi,chani,fi,:) = 10*log10( mean(abs(as(:, times2saveidx, LFP.cond  == condi)).^2,3) ./ basepow ); %
                itpcz(subno,condi,chani,fi,:) = size(as(:,:, LFP.cond == condi),3)* abs(mean(exp(1i*angle(as(:, times2saveidx, LFP.cond == condi))),3)).^2; % itpcz = ntrial*itpc²
            end % end condition loop
            
        end % end frequencies loop
        
    end% end channel loop
    
end% subject loop


%% Experimental condition coding

% 1 = 0_59 : 0 congruent
% 2 = 0_62 : 0 incongruent
% 3 = 1_59 : 1 congruent
% 4 = 1_62 : 1 incongruence
% 5 = 100_59 : 100 congruent
% 6 = 100_62 : 100 incongruent

%% Ploting


% Grand-average power
% 1st pair of contact

figure(1), clf

climdb  = [-1.5 1.5];
contourf(times2save,frex,squeeze(mean(mean(squeeze(tf(:,:,1,:,:)) ,2) ,1)),60,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),log10(max_frex),4),'yticklabel',round(logspace(log10(min_frex),log10(max_frex),4)*10)/10,'xlim', [-0.3 2])
title('Averaged time-frequency power', 'FontSize', 22)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
cbh = colorbar;
set(cbh, 'YTick', [-1.5, 0, 1.5])
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
cbh.Label.String = 'dB  from baseline'
cbh.FontSize = 24
ax = gca;
ax.FontSize = 24;
% Set size of the figure
x0=10;
y0=10;
width=700;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])
colormap jet


% If you want to look per condition
figure(2), clf

climdb  = [-2 2];
for condi=1:6
    subplot(3,2,condi)
    contourf(times2save,frex,squeeze(mean(squeeze(tf(:,condi,1,:,:)) ,1)),60,'linecolor','none')
    set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),log10(max_frex),6)*10)/10,'xlim', [-0.3 2])
    title({'Time-frequency power of condition' num2str(condi)})
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    colorbar
    ylabel(colorbar, 'dB  from baseline')
    colormap jet
end

% Grand-average ITPCz
% 1st pair of contact

figure(3),clf

climdb  = [0 3];
contourf(times2save,frex,squeeze(mean(mean(squeeze(itpcz(:,:,1,:,:)) ,2) ,1)),60,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),log10(max_frex),6)*10)/10,'xlim', [-0.3 2])
title('Averaged ITPCz')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'ITPCz')
colormap jet ;
cbh = colorbar;
set(cbh, 'YTick', [0, 1.5, 3])
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
cbh.Label.String = 'ITPCz'
cbh.FontSize = 24
ax = gca;
ax.FontSize = 24;
% Set size of the figure
x0=10;
y0=10;
width=700;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])
colormap jet

% If you want to look ITPCz per condition 

figure(4), clf
climdb  = [0 4];
for condi=1:6
    subplot(3,2,condi)
    contourf(times2save,frex,squeeze(mean(squeeze(itpcz(:,condi,1,:,:)) ,1)),60,'linecolor','none')
    set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),log10(max_frex),6)*10)/10,'xlim', [-0.3 2])
    title(condiname(condi))
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    colorbar
    ylabel(colorbar, 'ITPCz')
    colormap jet
end


%% Export power/ITPCz data from a specific time-frequency window

% Definition of TF windows

% ---- Change time_windows and freq_windows values according to the time-frequency window you are interested in --- %
% You need to use the exact values that are in times2save or frex

time_windows = [ 0.08 0.46 ]; 
freq_windows = [ 2.2890 8.2312 ]; 

% Intialize vectors for time and frequency indices
timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

% find the time indices corresponding to time and frequency windows
for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

% find the frequency indices corresponding to time and frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end


% from now, use the tf matrix if you want to extract power values, or the itpcz matrix if you want to extract ITPCz values
for condi=1:size(tf,2) % loop over conditions
    
        % pointer to stats file
        statsfilename=(['statistics_file_power_cog_theta', num2str(condi), '.txt']);

        fid=fopen(statsfilename,'w');
        
    for subno=0:size(tf,1) % loop over subjects
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1)) num2str(time_windows(ti,2)) '_f ' num2str(freq_windows(fi,1)) num2str(freq_windows(fi,2)) '\t' ]);
                   
                else
                    
                    % get data from large window
                    tempdat = squeeze(tf(subno,condi,1,freqidx(fi,1):freqidx(fi,2),timeidx(ti,1):timeidx(ti,2)));
                    % find max point (must first vectorize otherwise you get max for each line) to get the data around that subject's specific peak in the big TF window
                    [junk,maxpoint]=max(tempdat(:));
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    % get the indices from the larger matrix, not the indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    fprintf(fid,'%g\t',mean(mean(squeeze(tf(subno,condi,1,maxF-1:maxF+1,maxT-6:maxT+6)),1),2)); % 6 so as to have a proportional number of frequency (3) and timepoints in the subject TF window
                    
                end
               
            end
                        
        end
        
        fprintf(fid,'\n'); % careful this is critical for good organizations of lines!
    end % end subject loop
    fclose(fid);
end % end condition loop

%% End