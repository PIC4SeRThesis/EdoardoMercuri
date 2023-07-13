import os
import cv2 as cv
import numpy as np
import mediapipe as mp
import matplotlib.pyplot as plt
import scipy.io as sio

# Define in a list the videos names to elaborate
process=['Run1', 'Run2', 'Run3', 'Run4', 'Run5', 'Run6', 'Run7', 'Run8', 'Run9', 'Run10', 'Run12', 'Run13', 'Run14', 'Run15']

# Variables used to define the outputs of the elaboration
show_light_det = True
save_undist_vid = True
show_img_segm = True
save_pose_vid = True
show_feet_pos = True
save_feet_pos_plt = True

# Check presence of output folders
src_folder = 'Video/'
out_folder = 'Video_elaboration/'
if not os.path.exists(out_folder):
        os.makedirs(out_folder)
out_folder_plot='mat_files/plot_data/'
if not os.path.exists(out_folder_plot):
    os.makedirs(out_folder_plot)

# Load a .mat file containing the time in which to stop the video elaboration
load last_frames mat file to save time in video elaboration
mat_PI = 'mat_files/last_frames.mat'
container = sio.loadmat(mat_PI)
last_frames = dict.fromkeys(process)
for runs in process:
    last_frames[runs] = float(container['last_frames'][runs][0])

# Manual overwrite of the previous times, in case the data is not present
last_frames = {'Run1':0, 'Run2':0, 'Run3':0, 'Run4':0, 'Run5':0, 'Run6':0, 'Run7':0, 'Run8':0, 'Run9':0, 'Run10':0, 'Run11':0, 'Run12':0, 'Run13':0, 'Run14':0, 'Run15':0}                   

# LED turn on detection coordinates and area radius
x_light = 61
y_light = 818
radius_check = 10

# This dict is used to save some time on the start of the elaboration, specifing when to start in seconds the LED turn on detections
light_detection_delay = {'Run1':0, 'Run2':0, 'Run3':0, 'Run4':0, 'Run5':0, 'Run6':0, 'Run7':0, 'Run8':0, 'Run9':0, 'Run10':0, 'Run11':0, 'Run12':0, 'Run13':0, 'Run14':0, 'Run15':0} 

# If for some reasons the LED turn on has not been detected by the elaboration, but it is detectable by video inspection by a human, the value can be manually overwritten; put 0 to do the detection.
light_det_over = {'Run1':0, 'Run2':0, 'Run3':0, 'Run4':0, 'Run5':0, 'Run6':0, 'Run7':0, 'Run8':0, 'Run9':0, 'Run10':0, 'Run11':0, 'Run12':0, 'Run13':0, 'Run14':0, 'Run15':0}  

# Starting rotation angle of the image, will be checked interactively later
angle=0

# Variable that will be used later
meters = 40    

# Dimension of the image after the barrel distorsion correction
x_out = 3424     # This value has been found empirically for 40 meters
y_ou t= 200

# Number of sections in which the image will be divided to perform the pose detection
division_factor = 8

# Length in pixels of each section
x_block = int(x_out/division_factor)   

# Barrel distortion correction function
D = np.array([0., 0., 0., 0.])
feq = 1450

def img_undistort(frame):
    global save_undist_vid

    # This matrix is needed for the function, and it features the coordinates of the point that after the correction will be the center of the image
    K = np.array([[feq ,  0. , int(x_out/2 - x_center_offset)],
                  [ 0. , feq , int(y_out/2 - y_center_offset)],
                  [ 0. ,  0. , 1.  ]])
    Knew= K.copy()
    
    # Actual function responsible for the barrel distorion correction
    img_undist = cv.fisheye.undistortImage(frame, K, D, Knew=Knew)
    if save_undist_vid:
        out_undist.write(img_undist)            
    return img_undist

# Function used to get coordinates of mouse click on an image, used in LED detection area check
def light_onclick(event):
    global x_delta_input
    global y_delta_input
    global cid
    global fig
    global x_light
    global y_light
    if event.xdata != None and event.ydata != None:
        x_delta_input, y_delta_input = round(event.xdata), round(event.ydata)        
        x_light = x_light + (x_delta_input-radius_check)
        y_light = y_light + (y_delta_input-radius_check)
        fig.figure.canvas.mpl_disconnect(cid)
        plt.close()     
    return

# Function called when the image is closed, confirming the LED detection area
def light_onclose(event):
    global light_area_confirm
    light_area_confirm = True
    print("Light area detection confirmed")
    return

# LED detection area function
def light_area(frame):     
    global x_light
    global y_light
    global light_area_confirm
    global cid
    global fig
    light_area_confirm = False

    # Show LED detection area and call check confirmation function above
    while light_area_confirm == False:
        print("Click on the LED")
        area = frame[y_light-radius_check:y_light+radius_check, x_light-radius_check:x_light+radius_check]        
        fig = plt.imshow(area)
        cid = fig.figure.canvas.mpl_connect('button_press_event', light_onclick) 
        cid = fig.figure.canvas.mpl_connect('close_event', light_onclose)
        plt.suptitle(curr_run)
        plt.grid(True)
        plt.show()
    return light_area_confirm 

# LED turn-on detection function
def light_detection(frame): 
    global counter
    global light_detected
    global show_light_det
    global reference
    global last_frames
    # The first time the function gets called, the reference frame gets set
    if reference is None: 
        reference = frame[y_light-radius_check:y_light+radius_check, x_light-radius_check:x_light+radius_check]
        reference = cv.cvtColor(reference, cv.COLOR_BGR2GRAY)
        reference = cv.GaussianBlur(reference,(5,5),0)

    # The current frame getss cropped, color converted and gaussian blurred
    current = frame[y_light-radius_check:y_light+radius_check, x_light-radius_check:x_light+radius_check]
    current = cv.cvtColor(current, cv.COLOR_BGR2GRAY)
    current = cv.GaussianBlur(current,(5,5),0)
    
    # The color difference gets calculated
    diff=cv.absdiff(current,reference)    

    # A threshold on the difference highlights if there are major difference, this happens when the LED turn on
    _,thresholded=cv.threshold(diff,10,255,cv.THRESH_BINARY)

    # Enlarge the detected difference
    dilated=cv.dilate(thresholded,None,iterations=3)

    # Define the contours coordinates of the detected area
    contours, _ = cv.findContours(dilated, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)    

    # The contours are calculated only if the turn-on detection is successful, in this case, show the result, and return the number of the current frame
    if contours!=():
        if show_light_det:
            _, axs = plt.subplots(2, 2)
            axs[0, 0].imshow(reference[:,:])
            axs[0, 0].set_title("Input Image")
            axs[0, 1].imshow(diff[:,:])
            axs[0, 1].set_title("Difference")
            axs[1, 0].imshow(thresholded[:,:])
            axs[1, 0].set_title("Threshold applied")
            axs[1, 1].imshow(dilated[:,:])
            axs[1, 1].set_title("Dilatation")
            plt.tight_layout()
            plt.show()
        print("Light turn on detected at "+ str(float(f'{counter / fps:.3f}')))
        light_detected=counter        
        if last_frames[curr_run]==0:
            last_frames[curr_run]=(total_frames - light_detected) / fps
        return
    else:
        return

# Check rotation angle function, displays an image with grid to check is if the value specified alignes the image with the horizon
def rotation_angle(frame):
    global angle
    confirm=False
    while confirm==False:
        grid=rotate_image(frame, angle)
        print("Press any letter to confirm, or enter a new angle of rotation in degrees, positive is counterclockwise")
        plt.imshow(grid)
        plt.suptitle(curr_run + " " + str(angle) + "°")
        plt.grid(True)
        plt.show()
        try:
            rot_input=float(input())
            angle=rot_input
        except ValueError:
            print("Rotation confirmed")
            confirm=True
    return confirm

# Function to rotate the image
def rotate_image(image, angle):
    image_center = tuple(np.array(image.shape[1::-1]) / 2)
    rot_mat = cv.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv.warpAffine(image, rot_mat, image.shape[1::-1], flags=cv.INTER_LINEAR)
    return result

# Function used to get coordinates of mouse click on an image, used to find the center marker placed at 20 meters, that will be the center of the image after correction
def center_onclick(event):
    global cid
    global fig
    global x_center
    global y_center
    radius=100
    if event.xdata != None and event.ydata != None:
        x_center_input, y_center_input = round(event.xdata), round(event.ydata)  
        x_center=x_center+(x_center_input-radius)
        y_center=y_center+(y_center_input-radius-90)
        fig.figure.canvas.mpl_disconnect(cid)
        plt.close()     
    return

# Function called when the image is closed, confirming the center position coordinates
def center_onclose(event):
    global center_area_confirm
    center_area_confirm=True
    print("Center area detection confirmed")
    return

# Check center position for barrel distortion
def center_pos(frame):     
    global x_center
    global y_center
    global y_crop_upper_limit
    global y_crop_lower_limit
    global black_border_right
    global black_border_left
    global x_center_offset
    global y_center_offset
    global cid
    global fig
    global center_area_confirm
    radius=100
    center_area_confirm=False
    while center_area_confirm==False:
        area=frame[y_center-radius:y_center+radius, x_center-radius:x_center+radius]
        print("Press just below the center cone, that area must be at x=100, y=190 coordinates")
        fig=plt.imshow(area)
        cid=fig.figure.canvas.mpl_connect('button_press_event', center_onclick) 
        cid=fig.figure.canvas.mpl_connect('close_event', center_onclose)
        plt.suptitle(curr_run)
        plt.grid(True)
        plt.show()
        x_center_offset = x_center - x/2
        y_center_offset = y_center - y/2
        y_crop_upper_limit = int((y-y_out)/2 + y_center_offset)
        y_crop_lower_limit = int((y-y_out)/2 + y_center_offset+y_out)
        black_border_right = int((x_out-x)/2 + x_center_offset)
        black_border_left  = int((x_out-x)/2 - x_center_offset)
    return center_area_confirm 

# Pose detection function configuration
mp_drawing_styles = mp.solutions.drawing_styles
mp_drawing = mp.solutions.drawing_utils
mp_pose = mp.solutions.pose
pose_landmark = mp_pose.PoseLandmark
pose_config = mp_pose.Pose(static_image_mode=False, min_detection_confidence=0.75,min_tracking_confidence=0.2)

coord_mult=meters/x_out 

def pose_detection(frame):
    global save_pose_vid
    global index_progression
    img_block_RGB = cv.cvtColor(frame, cv.COLOR_BGR2RGB)

    # Actual pose detection function
    pose_info = pose_config.process(img_block_RGB)
    if pose_info.pose_landmarks:      

        # Change section of thewhere to perform the pose detection if one of the feet is near left border, in the last fifth on le left side 
        if (pose_info.pose_landmarks.landmark[31].x <= 0.2 or pose_info.pose_landmarks.landmark[32].x <= 0.2) and block_index>1:    

            # This variable determines which section is being processed to the pose detectino function, ranges from 1 to 15, which is 8 sections plus the 7 in between, which are composed by half of one section and half of the next one
            index_progression+=1

        if save_pose_vid:
        
            # Function ot draw the "stick-man" on top of the image
            mp_drawing.draw_landmarks(image = frame, 
                landmark_list = pose_info.pose_landmarks, connections = mp_pose.POSE_CONNECTIONS,
                landmark_drawing_spec = mp_drawing.DrawingSpec(color = (0,0,255), circle_radius = 1),
                connection_drawing_spec = mp_drawing.DrawingSpec(color = (49,125,237), circle_radius = 1))

            # Save the frame in a video
            out_pose.write(frame)

        # Get the position expressed in meters from 0 to 40, knowing the feet position in the section, size of each section, and number of the section
        x_LF_det = round(40 - ((pose_info.pose_landmarks.landmark[32].x * x_bloc k + (block_index-1) * x_block )* coord_mult), 2)
        x_RF_det=round(40 - ((pose_info.pose_landmarks.landmark[31].x * x_block + (block_index-1) * x_block) * coord_mult), 2)

    # If pose detection did not succeede, place the position to 0 and still save the frame without the "stick-man" in order to not skipt time frames on the video output which would desync the data from the video
    else:
        if save_pose_vid:
            out_pose.write(frame)
        x_LF_det=0
        x_RF_det=0
    return x_LF_det, x_RF_det

# Start of the big loop that will call all the functions above
for curr_run in process:

    # Load the Run# video file
    src = cv.VideoCapture(src_folder+curr_run+'.MP4')

    # Get video properties
    fps = src.get(cv.CAP_PROP_FPS)
    total_frames = src.get(cv.CAP_PROP_FRAME_COUNT)
    x = int(src.get(cv.CAP_PROP_FRAME_WIDTH )) # x pixel dimension
    y = int(src.get(cv.CAP_PROP_FRAME_HEIGHT)) # y pixel dimension
    
    x_center=int(x/2)
    y_center=int(y/2)

    fps_out_undist=fps
    fps_out_pose=fps
    
    # Image segmentation initialization
    index_progression = 0

    # LED turn-on detection initialization
    confirm_light_area = False
    reference = None
    light_detected= 0
    
    # Rotation initialization
    confirm_rotation_angle = False

    # Check center position for barrel distorsion initialization
    confirm_center_pos = False
    
    # Barrel distorsion video output initialization
    if save_undist_vid:
        fourcc = cv.VideoWriter_fourcc(*'mp4v')
        out_dist_corr=out_folder+curr_run+'_undistorted.mp4'
        out_undist = cv.VideoWriter(out_dist_corr,fourcc,fps_out_undist, (x_out,y_out))   

    # Pose detection video output initialization
    if save_pose_vid:
        fourcc = cv.VideoWriter_fourcc(*'mp4v')
        out_pose_det=out_folder+curr_run+'_pose.mp4'
        out_pose = cv.VideoWriter(out_pose_det,fourcc,fps_out_pose, (x_block,y_out))

    # Progress counter initialization
    counter = 0
    avanzamento = 0
    avanzamento_sec = 0

    # Output vectors initialization
    x_LF=[]
    x_RF=[]
    t=[]

    # Elaboration of each frame of the video
    while src.isOpened():    
        ret, frame = src.read()

        # Exit if the video file ended
        if not ret:
            print("Can't receive frame (stream end?). Exiting ...")
            break

        # Exit if the specified last frame is reached
        if light_detected !=0 and counter >= light_detected + last_frames[curr_run] * fps:
            print("Reached last step frame. Exiting ...")
            break

        # Counter increase
        counter += 1

        # LED turn-on
        if light_detected == 0:
            if light_det_over[curr_run] !=0 : # if overwrite active
                if counter >= light_det_over[curr_run] * fps:
                    light_detected = counter
            elif counter >= (light_detection_delay[curr_run]*fps): # LED turn-on detection starts after manual delay

                # Check LED area coordinates
                if confirm_light_area == False:
                    confirm_light_area=light_area(frame)
                else:
                    # Call LED turn-on detection function
                    light_detection(frame)

            # Progress feedback on the terminal
            q_sec = float(f'{counter/fps:.1f}')
            if q_sec > avanzamento_sec:
                avanzamento_sec = q_sec
                print(str(avanzamento_sec) + " seconds elapsed")

        # When turn-on is detected, elaborate the video
        else:
            
            # Image rotation angle
            if confirm_rotation_angle==False:
                confirm_rotation_angle=rotation_angle(frame)
            frame = rotate_image(frame, angle)

            # Image cropped to the middle vertically to reduce workload
            if confirm_center_pos==False:
                confirm_center_pos=center_pos(frame)
            frame = frame[ y_crop_upper_limit : y_crop_lower_limit , : ] 

            # Add black borders to the left and right to avoid cropping the image after barrel distortion correction
            frame = cv.copyMakeBorder(frame, 0, 0, black_border_left, black_border_right, cv.BORDER_CONSTANT, 0)

            # Barrel distortion correction
            frame = img_undistort(frame)

            # Frame segmentation preparing for pose detection
            block_index = division_factor - (index_progression * 0.5)  
            frame = frame[ : , int(x_block * (block_index-1)) : int(x_block * block_index)]

            # Pose detection function
            x_LF_det , x_RF_det = pose_detection(frame)

            # Position saved in the output vectors
            x_LF.append(x_LF_det)
            x_RF.append(x_RF_det)

            # Time saved in the output vector
            t.append(float((counter-light_detected)/fps))
        
            # Progress feedback on the terminal
            q=int(((counter-light_detected)*100) / (last_frames[curr_run]*fps))
            if q>avanzamento:
                avanzamento = q
                print(avanzamento,"%"," video elaborato")  

    if save_pose_vid:
        out_pose.release()
    if save_undist_vid:
        out_undist.release()

    # Save the output vectors as .mat files
    outdict={"x_RF":x_RF,"x_LF":x_LF,"t":t}
    sio.savemat(out_folder_plot+curr_run+'_txLxR.mat', outdict)

    # Show the positions to check if the pose detection succeeded
    if save_feet_pos_plt or show_feet_pos:
        fig, ax = plt.subplots()
        ax.scatter(t,x_LF,s=1,label="Left foot pos")
        ax.scatter(t,x_RF,s=1,label="Right foot pos")
        ax.set(xlabel='time (s)', ylabel='position (m)',title='Detected positions')    
        ax.legend()
        ax.grid()    
        plt.suptitle(curr_run)
        if save_feet_pos_plt:
            out_feet_pos_plt=out_folder+curr_run+'_feet_pos.png'
            plt.savefig(out_feet_pos_plt, bbox_inches='tight', dpi=1200)
        if show_feet_pos:
            plt.show()