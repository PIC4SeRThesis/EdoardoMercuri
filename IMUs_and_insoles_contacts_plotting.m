% Define a vector specifying which runs are meant to be processed
process=1:15;

% If true show the moments in which the foot is on the ground, 
% highlighted in red
show_contacts=true;

% Load insoles data struct
load('mat_files\010_results.mat')

% Get field names from insoles data struct
resultsEntries = fieldnames(results);

% Start a big for-loop to show a plot for each entry in 
% the "process" vector
for i=process

    % Inside the for-loop 2 figures are created, one with accelerations 
    % and one with angular velocities, with 6 plots each, 3 axes per 
    % foot. Only the first two plot for the acceleration figure and 
    % the first for the angular velocity are listed below.
    
    % Prepare a string with the shape "RunX" with X equal to the 
    % number of the current run
    Run=strcat('Run',int2str(i));
    
    % Load the corresponding IMUs log, it will load a struct with 
    % the name "inertia"
    load(strcat('mat_files\',Run,'.mat'))
    
    % Get the fields names in the IMUs struct: the nodes and 
    % gateway names
    ProMoveNames = fieldnames(inertia);
    
    % Store the data from rightsensor, i.e. "node_669"
    Data_R=inertia.(ProMoveNames{1}).data;  
    
    % Store the data from left sensor, i.e. "node_670"
    Data_L=inertia.(ProMoveNames{2}).data;  
    
    % Get the timestamp in which the port 1 got set to "high", 
    % the first index is the number of the registered event, 
    % the second index is the field of the timestamps, which 
    % is the first column.
    trigger_timest=inertia.(ProMoveNames{3}).data(1,1);
    
    % Insoles data is stored in such a way that all the data 
    % is referred to the port 1 set to "high" as origin of the 
    % timeline, so the timestamp registered by the gateway needs 
    % to subtracted the timestamps of the data of the IMUs so 
    % that the resulting 0 is th esame as the insoles data.
    % As before, the second index is the column of the 
    % timestamps, which is 1
    Timest_R=Data_R(:,1)-trigger_timest;    
    Timest_L=Data_L(:,1)-trigger_timest;
    
    % The IMUs data is sampled on 1 kHz of frequency, and needs 
    % to be downsampled to 100 Hz to emulate a cheaper sensor.
    Timest_R_DS = downsample(Timest_R,10);
    Data_R_DS   = downsample(Data_R,10);
    Timest_L_DS = downsample(Timest_L,10);
    Data_L_DS   = downsample(Data_L,10);   
    
    % Sign modifications
    Data_R_DS(:,2)=-Data_R_DS(:,2); % Vertical acc
    Data_L_DS(:,6)=-Data_L_DS(:,6); % Anteroposterior gyr
    Data_L_DS(:,7)=-Data_L_DS(:,7); % Mediolateral gyr
    
    % Store data of the insoles of the current, only if requested.
    % Every Run has a struct for the right foot (RF) and the 
    % left foot (LF).
    if show_contacts==true
        Data_PIR=results.(Run).RF;
        Data_PIL=results.(Run).LF;
    
        % Each foot structure has a field containing a vector with
        % all the initial contacts timestamps, and a separate for 
        % the final contact timestamps.
        % This data needs to be divided by the sampling frequency,
        % 200 Hz, to find the timings in seconds.
        Data_PIR.IC=Data_PIR.IC/200;
        Data_PIR.FC=Data_PIR.FC/200;
        Data_PIL.IC=Data_PIL.IC/200;
        Data_PIL.FC=Data_PIL.FC/200;
    end
    
    % Acceleration figure
    figure('Name',append(Run,' Acceleration 100 Hz'))
    
    %% Vertical Acceleration (positive upwards)
    subplot(2,3,1);
    
    % The 2 on the second index of the data is taking the x axis
    % of the accelerometer
    plot(Timest_L_DS, Data_L_DS(:,2)) 
    xlim([0 Timest_L_DS(end)])
    
    % If requested show every contacts as an area, combining the 
    % initial and final contact vectors.
    if show_contacts==true
        for j=1:1:length(Data_PIL.IC)
            hold on
            area([Data_PIL.IC(j) Data_PIL.FC(j)],[160 160], ...
            'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
            hold on
            area([Data_PIL.IC(j) Data_PIL.FC(j)],[-160 -160] ...
            ,'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
        end
        hold off
        xlim([Data_PIL.IC(1) Data_PIL.FC(end)])
    end
    grid on
    title('Vertical Left')
    ylabel('m/s^2')
    xlabel('s')
    ylim([-160 160])
    
    
    subplot(2,3,4);
    % On the right foot it can be noted that the data is taken 
    % with opposite sign, to have positive values when upwards.
    plot(Timest_R_DS, -Data_R_DS(:,2)) 
    xlim([0 Timest_R_DS(end)])    
    if show_contacts==true
        for j=1:1:length(Data_PIR.IC)
            hold on
            area([Data_PIR.IC(j) Data_PIR.FC(j)],[160 160], ...
            'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
            hold on
            area([Data_PIR.IC(j) Data_PIR.FC(j)],[-160 -160], ...
            'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
        end    
        hold off
        xlim([Data_PIR.IC(1) Data_PIR.FC(end)])
    end
    grid on
    title('Vertical Right')
    ylabel('m/s^2')
    xlabel('s')
    ylim([-160 160])
    
    % [...]
    
    %% Gyroscope figure
    figure('Name',append(Run,' Angular velocity 100 Hz'))
    
    %% Vertical axis Angular velocity
    % Orientation of vectors left as they are, so a clockwise
    % rotation of the right foot is the same sign of a 
    % counterclockwise rotation on the left foot.
    subplot(2,3,1);
    plot(Timest_L_DS, Data_L_DS(:,5))
    xlim([0 Timest_L_DS(end)])
    if show_contacts==true
        for j=1:1:length(Data_PIL.IC)
            hold on
            area([Data_PIL.IC(j) Data_PIL.FC(j)],[1600 1600], ...
            'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
            hold on
            area([Data_PIL.IC(j) Data_PIL.FC(j)],[-1600 -1600], ...
            'FaceColor','#D95319','EdgeColor','none', ...
            'FaceAlpha',0.3)
        end    
        hold off
    xlim([Data_PIL.IC(1) Data_PIL.FC(end)])
    end
    grid on
    title('Vertical Left')
    ylabel('rad/s')
    xlabel('s')
    ylim([-1600 1600])
    % scrollplot;
end