% This code illustrates the SVD compression of the signal wave of ECG
% signal
% Authors: Amol Shah, Mantri Krishna Sri Ipsit and Shlok Vaibhav Singh
% (Group 1) of EE338 Autumn 2020 
clc;
clear all;


%%%%%%%%%%%%% Parameters that can be modified
no_of_cycles = 1000;   % No of R-to-R peak cycles 
q = 3;    % Increase q for better reconstruction but lesser compression

%%%%%%%%
% Can modify the dataset by changing 103 to desired daatset number of
% MIT-BIH, check net for valid values
[sig, Fs, tm] = rdsamp('mitdb/103', 1);
[record, ann] = rdann('mitdb/103','atr'); % R peaks data



Rpeak_samples =  record(find(ann=='N'));
no_of_cycles = min(no_of_cycles+1,length(Rpeak_samples))-1;
Rpeak_samples = Rpeak_samples(1:no_of_cycles+1);
cycle_durations = record(3:no_of_cycles+2) - record(2:no_of_cycles+1);
avg_duration = ceil( mean(cycle_durations));

%Period normalisation transformation
A = zeros(no_of_cycles,avg_duration);
for i = 1:no_of_cycles
    old_values = sig(Rpeak_samples(i):Rpeak_samples(i+1));
    cycle_duration = cycle_durations(i);
    for j=1:avg_duration
        rj = (j-1)*(cycle_duration-1)/(avg_duration-1)+1;
        
        jstar = floor( rj);
        A(i,j) = old_values(jstar) + (old_values(jstar+1) - old_values(jstar))* (rj-jstar);
    end
end

% SVD Extraction
[U,S,V] = svd(A);
extracted_U = U(:,1:q);
extracted_V = V(:,1:q);
extracted_S = diag(S);
extracted_S = extracted_S(1:q);

% Reconstruction 
Ahat = zeros(no_of_cycles,avg_duration);
for i = 1:q
    temp = extracted_U(:,i)*extracted_S(i)*transpose(extracted_V(:,i));
    Ahat = Ahat + temp;
end

reconst_sig= [];
% Note that we used only Ahat and period duration information for
% reconstruction 
for i=1:no_of_cycles
    
    old_values = zeros(avg_duration +1,1);
    old_values(1:end-1) = Ahat(i,:);
    old_values(end) = Ahat(i,end-1);
    new_val = zeros(cycle_durations(i),1);
    
    for j=1:cycle_durations(i)
        rj = (j-1)*(avg_duration-1)/(cycle_durations(i)-1)+1;
        
        jstar = floor(rj);
        new_val(j) = old_values(jstar) + (old_values(jstar+1) - old_values(jstar))* (rj-jstar);
    end
    reconst_sig = cat(2,reconst_sig,transpose(new_val));
end

% Comparision of original and recontructed
original_sig = transpose(sig(Rpeak_samples(1):Rpeak_samples(end)-1));
duration = 10;  % Plot duration shown
plot(tm(1:360*duration),original_sig(1:360*duration));
hold on;
plot(tm(1:360*duration),reconst_sig(1:360*duration));
xlabel("Time")
ylabel("ECG Signal in mV")

% PRD calculation
temp = sum((original_sig - reconst_sig).^2);
temp2 = sum(original_sig.^2);
PRD = sqrt(temp/temp2) *100



