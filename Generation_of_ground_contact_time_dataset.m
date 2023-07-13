% Prepare 2 empty vectors, populated at the end of each run 
% elaboration with new samples
y=[];
x=[];

% Moving average filter periodicity
N=6;

% Length of the samples, given that the data is sampled at 
% 100 Hz, the time windows is 500 ms 
T=50;

% Positive threshold to filter the local maxima
threshold=400;

% Define a vector specifying which runs are meant to be processed
process=1:8;

% Load insoles data struct
load('mat_files\001_results.mat')

% Get field names from insoles data struct
resultsEntries = fieldnames(results);

% Start a big for-loop to show augment the "x" and "y" matrices 
% for each entry in the "process" vector
for i=process

    % The following code is the for-loop, when the functions employed 
    % are similar if not equal to the previously detailed code they 
    % will not be commented.

    % Load insoles data
    Run=append('Run',int2str(i));
    
    Data_PIR=results.(Run).RF;
    Data_PIL=results.(Run).LF;
    
    Data_PIR.IC         = Data_PIR.IC/200;
    Data_PIR.FC         = Data_PIR.FC/200;
    Data_PIL.IC         = Data_PIL.IC/200;
    Data_PIL.FC         = Data_PIL.FC/200;
    
    % Add two new fields to the structure: duration and middle 
    % timestamp of each contact
    Data_PIR.DurationC  = Data_PIR.FC-Data_PIR.IC;
    Data_PIL.DurationC  = Data_PIL.FC-Data_PIL.IC;
    Data_PIR.TimeC      = (Data_PIR.FC+Data_PIR.IC)/2;
    Data_PIL.TimeC      = (Data_PIL.FC+Data_PIL.IC)/2;
    
    % Preallocate memory to a matrix used to store both left and 
    % right data in a compact way
    steps_temp = zeros(length(Data_PIR.IC)+length(Data_PIL.IC),2);
    
    % Store middle timestamps in the first column, and contact time 
    % in the second
    steps_temp(:,:) = [Data_PIR.TimeC' Data_PIR.DurationC';
                       Data_PIL.TimeC' Data_PIL.DurationC'];
    
    % Sort the matrix using the first column in ascending order, 
    % save also the sorting order
    [steps_temp(:,:), R_L_sort_y] = sortrows(steps_temp(:,:),1);
    
    
    load(strcat('mat_files\',Run,'.mat'))
    
    ProMoveNames = fieldnames(inertia);
    Data_R=inertia.(ProMoveNames{1}).data;
    Data_L=inertia.(ProMoveNames{2}).data;
    
    trigger_timest=inertia.(ProMoveNames{3}).data(1,1);
    
    % In addition to reset the timestamps here the values 
    % before zero get discarded
    [~,cutR]=min(abs(Data_R(:,1)-trigger_timest));
    [~,cutL]=min(abs(Data_L(:,1)-trigger_timest));
    Data_R=Data_R(cutR:end,:);
    Data_L=Data_L(cutL:end,:);
    Data_R(:,1)=Data_R(:,1)-trigger_timest;
    Data_L(:,1)=Data_L(:,1)-trigger_timest;    
    Timest_R=Data_R(:,1);
    Timest_L=Data_L(:,1);
    
    % The two matrices are never the same size because of 
    % start and stop mismatch, but here it is detected if 
    % one of the twwo is significally less than the other 
    % (ratio bigger than 1.5).
    if length(Timest_R)/length(Timest_L)>1.5
        fprintf(strcat("Bad Left data in Run ", ...
        int2str(i)),", filling missing samples...")
        
        % This can happen on the final samples where it is 
        % not crucial to have data since insoles data is not 
        % present, but in order to not occur in problems 
        % during further elaborations, here the lost samples
        % are filled with the average of the two nearest 
        % values.
        Data_L_temp=zeros(length(Timest_R),13);
        for j=1:length(Timest_R)
            if ismember(Timest_R(j),Timest_L)
                [~,id]=min(abs(Timest_L-Timest_R(j)));
                Data_L_temp(j,:)=Data_L(id,:);
            else
                [~,id]=min(abs(Timest_L-Timest_R(j)));
                if id<length(Data_L(:,1))
                    Data_L_temp(j,:)=[Timest_R(j) ...
                    (Data_L(id,2:end)+Data_L(id+1,2:end))/2];
                else
                    Data_L_temp(j,:)=[Timest_R(j) Data_L(id,2:end)];
                end
            end
        end
        Data_L=Data_L_temp;
        Timest_L=Data_L(:,1);
    
        % Same if the right data has lost some samples.
    elseif length(Timest_L)/length(Timest_R)>1.5
        % [...]
    end
    
    Timest_R_DS = downsample(Timest_R,10);
    Timest_L_DS = downsample(Timest_L,10);
    Data_R_DS   = downsample(Data_R,10);
    Data_L_DS   = downsample(Data_L,10);
    
    Data_R_DS(:,2)=-Data_R_DS(:,2); % Vertical acc
    Data_L_DS(:,6)=-Data_L_DS(:,6); % Anteroposterior gyr
    Data_L_DS(:,7)=-Data_L_DS(:,7); % Mediolateral gyr
    
    % Save gyroscope mediolateral axis for local maxima 
    % identification
    Data_R_gyrmed = zeros(length(Data_R_DS(:,7)),1);
    Data_L_gyrmed = zeros(length(Data_L_DS(:,7)),1);
    
    % Moving average filter, N periodicity, right foot
    for k = 1+N:1:length(Data_R_DS(:,7))-N
        for l = -N:1:N
            Data_R_gyrmed(k)=Data_R_gyrmed(k) + ...
            Data_R_DS(k+l,7);
        end
        Data_R_gyrmed(k)=Data_R_gyrmed(k)/(2*N+1);
    end
    
    % Left foot
    for k = 1+N:1:length(Data_L_DS(:,7))-N
        % [...]
    end
    
    % Right foot mediolateral angular velocity local maxima, 
    % filtered by a minimum a positive "threshold".
    [Data_R_peaks,Data_R_peaks_ids] = findpeaks(Data_R_gyrmed);
    l=0;
    for k = 1:1:length(Data_R_peaks)
        if Data_R_peaks(k-l)<threshold
            % Empy the value if it is below the threshold.
            Data_R_peaks(k-l)=[];
            Data_R_peaks_ids(k-l)=[];
            l=l+1;
        end
    end    
    
    % Find timestamps of the local maxima, discard them if they
    % are before first step's initial contact time.
    for k = 1:1:length(Data_R_peaks)
        if Timest_R_DS(Data_R_peaks_ids(k))>Data_PIR.IC(1)
            Data_R_peaks = Data_R_peaks(k-1:end);
            Data_R_peaks_ids = Data_R_peaks_ids(k-1:end);
            break
        end
    end
    
    % Discard the local maxima after the last step's final 
    % contact time.
    for k = 1:1:length(Data_R_peaks)
        if Timest_R_DS(Data_R_peaks_ids(k))>Data_PIR.FC(end)
            Data_R_peaks=Data_R_peaks(1:k);
            Data_R_peaks_ids=Data_R_peaks_ids(1:k);
            break
        end
    end
    
    % Find time ids in the middle of the local maxima
    Data_R_center_time_ids=zeros(length(Data_R_peaks)-1,1);
    for k = 1:1:length(Data_R_center_time_ids)
        Data_R_center_time_ids(k) = ...
        round((Data_R_peaks_ids(k+1) + Data_R_peaks_ids(k))/2);
    end
    
    % Save the relative timestamps
    Data_R_center_timest = Data_R_DS(Data_R_center_time_ids);
    
    % Left foot data elaboration
    % [...]
    Data_L_center_timest = Data_L_DS(Data_L_center_time_ids);
    
    % Sorting the timestamps should lead to the same sorting 
    % order of the steps timings.
    [~, R_L_sort_x] = sortrows([Data_R_center_timest; ...
    Data_L_center_timest],1);
    
    if R_L_sort_x~=R_L_sort_y
        fprintf("Steps sorting mismatch between x and y.")
    end
    
    % Right foot plot.
    figure(i)
    % Show processed gyroscope mediolateral axis with its 
    % local maxima.
    plot(Timest_R_DS, Data_R_gyrmed , ... 
    Timest_R_DS(Data_R_peaks_ids), Data_R_peaks,'r*')
    hold on
    
    % Show middle points of the local maxima.
    scatter(Data_R_center_timest,1000,"b*")
    
    % For each step show an area when there is a contact.
    for l=1:1:length(Data_PIR.IC)
        hold on
        area([Data_PIR.IC(l) Data_PIR.FC(l)],[1000 1000] ...
        ,'FaceColor','#D95319','EdgeColor','none', ...
        'FaceAlpha',0.3)
        hold on
        area([Data_PIR.IC(l) Data_PIR.FC(l)],[-1000 -1000] ...
        ,'FaceColor','#D95319','EdgeColor','none', ...
        'FaceAlpha',0.3)
    end    
    title("Right foot")
    xlim([Data_PIR.IC(1)-2 Data_PIR.FC(end)+2])
    grid on
    
    % Allocate memory to a matrix where the sensors data 
    % will be stored.
    Data_R_windows=zeros(length(Data_R_center_time_ids),T+1,7);
    
    % For each middle point of the local maxima.
    for k = 1:1:length(Data_R_center_time_ids)
    
        % In each element of the first index, save a T+1x7 
        % matrix containing the data in the neighbourhood of 
        % the middle point of the local maxima. 
        % The first 7 columns contain timestamps, 3 
        % accelerometer axes and 3 gyroscope axes.
        Data_R_windows(k,:,:) = ...
        Data_R_DS(Data_R_center_time_ids(k)-T/2 : ...
        Data_R_center_time_ids(k)+T/2, 1:7);
    
        % Show just the saved gyroscope mediolateral axis to 
        % see if it is overlapping with the processed vector.
        hold on
        plot(Data_R_windows(k,:,1), Data_R_windows(k,:,7))
    end
    
    x_steps=Data_R_windows;
    
    % Left foot plot.
    figure(i+20)
    % [...]
    
    % Add the left foot's samples to the right foot's samples.
    x_steps=[ x_steps ; Data_L_windows];
    
    % Sort the first index with the sorting order of the 
    % timestamps.
    x_temp=x_steps(R_L_sort_x,:,:);
    
    % Save the step's contact times.
    y_temp=steps_temp(:,2);
    
    % Add the samples and labels to the global matrices, which 
    % are augmented each cycle of the for loop.
    x=[x;x_temp]; 
    y=[y;y_temp];   
end    

% Save the two matrices in two different files in another folder.
cd('..')
% save('Dataset_merge\1_y_CT.mat','y','-v7','-nocompression');
% save('Dataset_merge\1_x_CT.mat','x','-v7','-nocompression');