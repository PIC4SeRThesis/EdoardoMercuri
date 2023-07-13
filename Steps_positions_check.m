% Define a vector specifying which runs are meant to be processed
process = 1:15;

% Define a vector specifying which runs are meant to be checked manually, this can be useful if some of the runs were already checked in a previus check with this script
check_process = 1:15;

% Manually specify the number of the last video run
runs_number = 15;

% Mode function parameters
average_range=20;
mode_toll=0.02;     % meters

% Half width in pixel of the image shown to correct the foot position
check_w=200;

% Load the insoles results if present, this helps speed up the checking process
if isfile('mat_files\001_results.mat')
    load('mat_files\001_results.mat')
    resultsEntries = fieldnames(results);
else 
    results={};    
end

% Specify if and how much the PI is late (+) or early (-), due to problems from light detection in seconds for each run 
t_adjust_R=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
t_adjust_L=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% Vectors initialization
max_steps_per_side=18;
steps_time_R = zeros(runs_number, max_steps_per_side);
steps_time_L = zeros(runs_number, max_steps_per_side);
steps_temp = zeros(max_steps_per_side * 2, 3, runs_number);
steps_time = zeros(runs_number, max_steps_per_side * 2);
steps_pos  = zeros(runs_number, max_steps_per_side * 2);
steps_side = zeros(runs_number, max_steps_per_side * 2); % R=0, L=1

% Load a previous check, if present, will be reprocessed according to the specified runs in the first vectors
if isfile('mat_files\steps.mat')
    load('mat_files\steps.mat')
end

% Start of the for loop for each trial
for i = process
    
    Run=strcat('Run',int2str(i));
    
    % Load undistorted video
    v = VideoReader(strcat('Video_elaboration/',Run,'_undistorted.mp4'));
    fps = v.FrameRate;
    width = v.Width;

    % Load steps positions and time from pose detection, variables loaded: t, x_LF, x_RF
    load(strcat('mat_files\plot_data\',Run,'_txLxR.mat')); 

    % Load pressure inertial insoles data if present
    if isfield(results,(Run))
        % Preprocess data of the insoles
        Data_PIR=results.(Run).RF;
        Data_PIL=results.(Run).LF;
        
        % Sampling frequency 
        fs=results.(Run).fs;                

        # Load IMUs data
        load(strcat('mat_files\',Run,'.mat'));
        
        % Timestamps in milliseconds
        Data_PIR.IC         = Data_PIR.IC/fs-t_adjust_R(i);    
        Data_PIR.FC         = Data_PIR.FC/fs-t_adjust_R(i);
        Data_PIL.IC         = Data_PIL.IC/fs-t_adjust_L(i);
        Data_PIL.FC         = Data_PIL.FC/fs-t_adjust_L(i);
        
        % Contact durations
        Data_PIR.DurationC  = Data_PIR.FC-Data_PIR.IC;
        Data_PIL.DurationC  = Data_PIL.FC-Data_PIL.IC;
        
        % Middle timestamps of contacts
        Data_PIR.TimeC      = (Data_PIR.FC+Data_PIR.IC)/2;
        Data_PIL.TimeC      = (Data_PIL.FC+Data_PIL.IC)/2;

        % Assign the average times of contact to steps_time vectors
        steps_time_R(i,1:length(Data_PIR.TimeC)) = Data_PIR.TimeC;
        steps_time_L(i,1:length(Data_PIL.TimeC)) = Data_PIL.TimeC;

        % Build a temporary matrix assigning a 0 or a 1 in the third column, if the step is of the right foot or the left respectively
        steps_temp(:,:,i)=[steps_time_R(i,:)' zeros(max_steps_per_side,1) zeros(max_steps_per_side,1);
                           steps_time_L(i,:)' zeros(max_steps_per_side,1) ones(max_steps_per_side,1)];

        % Sort the rows of the matrix according to the first column of the average times of contact
        steps_temp(:,:,i)=sortrows(steps_temp(:,:,i));
        
        % Find middle timestamps of contacts in the pose detection time vector
        for j = 1:1:length(steps_temp(:,1,i))
            
            % Find the closest timestamp to the middle timestamp of j contact, save its index
            [~,id]=min(abs(t-steps_temp(j,1,i)));
            if id~=1
            
                % Collection of positions of the neighbourhood of the middle timestamp in a temporary vector, both feet are collected because the pose detection swaps them arbitarly around turnover
                step_temp = [x_RF(id-average_range : min(id+average_range , length(x_RF))) x_LF(id-average_range : min(id+average_range , length(x_LF)))]';
                
                % Custom mode function, 2 cm tollerance, average of the first 3 values
                steps_temp(j,2,i)=mode(step_temp,mode_toll,3);
            end
        end
        
        % Filtering out positions equal to 0 until the first non-zero value; successive zeros are holes to fix later
        j=1;
        while steps_temp(j,2,i) == 0
            j=j+1;
        end

        % The first non-zero value gets shifted in first position
        steps_temp(1:end-(j-1),:,i) = steps_temp(j:end,:,i);
        
        % Clean the leftovers from the shifting
        steps_temp(end-(j-2):end,:,i) = 0;

        % Add missing positions from: end-(j-1)+k
        c=false;
        k=1;
        while c==false
            fh = figure(i);
            fh.WindowState = 'maximized';

            % Show the positions of the pose detection
            scatter(t, x_RF, 4, "filled", "color", "#D95319")
            hold on
            scatter(t, x_LF, 4, "filled", "color", "#0072BD")
            hold on

            % Show the currently detected steps
            scatter(steps_temp(:,1,i), steps_temp(:,2,i), 400, "+", "color", "b")
            hold on

            % Linear regression of the steps overlayed to guide the inclusion of new steps
            X = [ones(length(steps_temp(1:end-(j-1),1,i)),1) steps_temp(1:end-(j-1),1,i)];
            b = X\steps_temp(1:end-(j-1),2,i);
            yCalc = b(1)+b(2)*t;
            plot(t,yCalc,"color","k")
            hold on

            % Labelled as 0 or 1 according to which foot it is
            labels = num2cell(steps_temp(:,3,i));
            text(steps_temp(:,1,i), steps_temp(:,2,i), labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')    
            ylim([-1 41])
            %  xlim([steps_temp(1,1,i)-3 max(steps_temp(:,1,i))+3])
            title(strcat(Run," manually add missing steps"))
            grid minor

            % Wait for mouse click on the graph to add missing steps, close the image to move to the next trial
            try [st,sp]=ginput(1);
            
                % Save the time and position of the new steps
                steps_temp(end-(j-1)+k,1,i)=st; 
                steps_temp(end-(j-1)+k,2,i)=sp;
                k=k+1;
                hold off
            catch
                c = true;
            end
        end

        % Sort the matrix to accomodate the new steps
        steps_temp(:,:,i)=sortrows(steps_temp(:,:,i));

        % Find the first 0 to know the vector size
        j=1;
        while steps_temp(j,2,i) == 0
            j=j+1;
        end        

        % Save times, positions and foot side to the steps_ vectors
        steps_time(i, 1:end-(j-1)) = steps_temp(j:end, 1, i)';
        steps_pos(i, 1:end-(j-1)) = steps_temp(j:end, 2, i)';
        steps_side(i, 1:end-(j-1)) = steps_temp(j:end, 3, i)';

    % If insoles data is not present, full manual mode
    else
        % Show time and position from pose detection,
        c=false;
        j=1;
        while c==false
            fh = figure(i);
            fh.WindowState = 'maximized';
    
            scatter(t,x_RF,4,"filled","color","#D95319")
            hold on
            scatter(t,x_LF,4,"filled","color","#0072BD")
            hold on
            scatter(steps_time(i,:),steps_pos(i,:), 400,"+","color","#7E2F8E")
            hold on
            
            if j==1 
                X = [ones(length(steps_time(i,1:j)),1) steps_time(i,1:j)'];
                b = X\steps_pos(i,1:j)';
            else
                X = [ones(length(steps_time(i,1:j-1)),1) steps_time(i,1:j-1)'];
                b = X\steps_pos(i,1:j-1)';
            end
            yCalc = b(1)+b(2)*t;
            plot(t,yCalc,"color","k")
            hold on
            ylim([-1 41])
            % xlim([1 12])
            title(strcat(Run," manually select steps"))
            grid minor
            
            % Manually select steps, close to confirm
            try [st,sp]=ginput(1);
                steps_time(i,j)=st;
                steps_pos(i,j)=sp;
                j=j+1;
                hold off
            catch
                c = true;
            end
        end        
    end
    
    % Foot side assignment
    step_pos_px = check_w;                 

    % Get video frame with the first foot on the ground
    Step_frame = read(v,round(steps_time(i,1) * fps));   

    % Draw a white vertical line on the calculated position
    Step_frame(:, width-step_pos_px, :) = 255;

    % Before displaying the frame, fill with black (value=0) pixels if the position is around the borders
    if width-step_pos_px <= check_w
        frame_fill = zeros(200, check_w - (width - step_pos_px), 3);
        Step_frame = [frame_fill Step_frame(:, 1:width-(step_pos_px-check_w), :)];
    elseif step_pos_px <= check_w
        frame_fill = zeros(200, check_w - step_pos_px, 3);
        Step_frame = [Step_frame(:,width-(step_pos_px+check_w):width,:) frame_fill];
    else

        % If not near a border, select a portion of the undistorted frame
        Step_frame = Step_frame(:, width - (step_pos_px + check_w):width - (step_pos_px - check_w), :);
    end

    % Show the frame
    f=figure(i);
    imshow(Step_frame)
    f.WindowState = 'maximized';
    title(strcat(Run," side selection"))
    
    % Expected click on the image right side for right foot and viceversa
    [x,~]=ginput(1);

    % The side vector is full of "0", assign the "1" for the left foot to the even or odd indexes according to the click on the image
    if x>check_w
        steps_side(i,2:2:end) = 1;
    else
        steps_side(i,1:2:end) = 1;
    end
    close(f)

    % Show final plot of steps and side assigned
    fh = figure(i+20);
    fh.WindowState = 'maximized';
    scatter(t,x_RF,4,"filled","color","#D95319")
    hold on
    scatter(t,x_LF,4,"filled","color","#0072BD")
    hold on
    scatter(steps_time(i,:),steps_pos(i,:),400,"+","color","k")
    labels = num2cell(steps_side(i,:));
    text(steps_time(i,:),steps_pos(i,:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')    
    ylim([-1 41])
    xlim([steps_time(i,1)-3 max(steps_time(i,:))+3])
    title(strcat(Run,' process end'))
    grid minor
    save("mat_files\steps.mat","steps_time","steps_pos","steps_side")
end

% After the missing steps are selected, their position needs to be checked
for i = check_process

    Run=strcat('Run',int2str(i));
    
    % Load undistorted video
    v = VideoReader(append(append('Video_elaboration/',Run),'_undistorted.mp4'));

    fps = v.FrameRate;
    width = v.Width;

    load(append('mat_files\plot_data\',append(Run,'_txLxR.mat'))); 
    d=false;
    
    % Process each step
    for j=1:length(steps_pos(i,:))

        % End of the process when a position = 0 is encountered
        if d==true
            steps_pos(i,j:end)=0;
            steps_time(i,j:end)=0;
            break
        end
        
        % End of the process when a time = 0 is encountered
        if steps_time(i,j)==0
            break
        end

        % Avoid problems with time and position limits
        if round((steps_pos(i,j)/40)*width)>width
            steps_pos(i,j)=39.95;
        end
        if round(steps_time(i,j)*fps)>v.NumFrames
            steps_time(i,j)=v.NumFrames/fps;
        end

        c = false;
        while c == false
        
            % Convert position from meters back to pixels
            step_pos_px = round((steps_pos(i,j)/40)*width);

            % Get video frame with a foot on the ground
            Step_frame = read(v,round(steps_time(i,j)*fps));

             % Draw a white vertical line on the calculated position
            Step_frame(:,width-step_pos_px,:) = 255;           
            
            % Show an image of the athlete. Before displaying the frame, fill with black (value=0) pixels if the position is around the borders
            if width-step_pos_px <= check_w
                frame_fill = zeros(200,check_w - (width - step_pos_px), 3);
                Step_frame = [frame_fill Step_frame(:, 1:width-(step_pos_px-check_w), :)];
            elseif step_pos_px <= check_w
                frame_fill = zeros(200,check_w - step_pos_px, 3);
                Step_frame = [Step_frame(:, width-(step_pos_px + check_w):width, :) frame_fill];
            else
                Step_frame = Step_frame(:, width - (step_pos_px + check_w):width - (step_pos_px -check_w), :);
            end
            f=figure(j);
            imshow(Step_frame)
            f.WindowState = 'maximized';
            title(strcat(Run,", step ",int2str(j)," position check"))
            % Expected click on the image to correct the detected position, close the window to confirm.
            try [x,~]=ginput(1);
                steps_pos(i,j)=((step_pos_px-(x-check_w))/width)*40; % Save the new position in meters
                if steps_pos(i,j)<0 || steps_pos(i,j)>40
                    steps_pos(i,j)=0;
                    d= true;
                end
                close(f)
            catch
                c = true;
            end
        end
    end

    % Show final plot of steps and side assigned with positions corrected
    fh = figure(i+20);
    fh.WindowState = 'maximized';
    scatter(t, x_RF, 4, "filled", "color", "#D95319")
    hold on
    scatter(t, x_LF, 4, "filled", "color", "#0072BD")
    hold on
    scatter(steps_time(i,:), steps_pos(i,:), 400, "+", "color", "k")
    labels = num2cell(steps_side(i,:));
    text(steps_time(i,:), steps_pos(i,:), labels, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')    
    ylim([-1 41])
    xlim([steps_time(i,1)-3 max(steps_time(i,:))+3])
    title(strcat(Run,' process end'))
    grid minor

    % Save the vectors each trial, to prevent lost work in case some problems occur in future trials
    save("mat_files\steps.mat","steps_time","steps_pos","steps_side")
end

% Custom mode funcition
function m = mode(lst, mode_toll, average_n)
    freq = zeros(length(lst),1);
    table = [lst freq];
    for i = 1:1:length(lst)
        for j = 1:1:length(lst)
            if abs(lst(i)-lst(j)) <= mode_toll
                table(i,2) = table(i,2) + 1;
            end
        end
    end
    table = sortrows(table,2);
    m = mean(table(end-average_n:end,1));
end