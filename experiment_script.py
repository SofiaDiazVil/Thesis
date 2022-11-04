from psychopy import visual, event, core, gui, constants, clock
import time, numpy, os, pandas
##start measuring time from the moment we start running the code 
timer = core.Clock()

## creat directory for Data
# set the directory
my_directory = os.getcwd()

# creat the folders if they dont exist already 
directory_to_write_to = my_directory + '/Data_ExpGent'
if not os.path.isdir(directory_to_write_to): #check if the file already exists 
    os.mkdir(directory_to_write_to) #creat a directory named data
    
## dialog box to collect the participants data 
# create a dictionary for the dialogue box 
InfoParticipant = {'Participant number':'1', 'Gender':["Male", "Female", "X"], 'Age':'23'}
   
## Present the dialogue box and check if the data file already exists 
# loop to check if a file already exists and if it does keep asking for a new file name (different participant number)
already_exists = True
while already_exists:
    # show the dialogue box
    infoDlg = gui.DlgFromDict(dictionary=InfoParticipant, title='Experiment Synchronization')
    
    # save the participants information (used later)
    participants_number = int(InfoParticipant['Participant number']) 
    gender = InfoParticipant['Gender']
    age = int(InfoParticipant['Age'])

    # creat the file name 
    file_name = directory_to_write_to + '/ecg_subj' + str(InfoParticipant['Participant number']) 
    
    # loop to see if a new file can be created or not    
    if not os.path.isfile(file_name): # check if file name exists 
        already_exists = False # name didn't exist so a new file can be made, the loop can stop 
    else: # name does exists so present a dialogue box informing this and ask for a new participants number 
        infoDlg = gui.Dlg(title = 'Error')
        infoDlg.addText('Try another participant number')
        infoDlg.show()


## constants 
# elements 
speedy = False #speed mode for code testing 
win= visual.Window(units= "norm", fullscr = True) #(size=[1000,700], units= "norm") # inititate window 
MessageOnSCreen = visual.TextStim(win, text = "OK") # text component for general text 
clock = core.Clock()
video = visual.MovieStim(win, filename = "") # video component for trials
example_rating = visual.ImageStim(win, image = "manekin.jpg", pos = (0.0, -0.4), size = (1,1.2)) # rating example for instructions
N_Trials = 2 # number of videos is 2, one erotic one neutral
trials = numpy.ones((N_Trials,9)) * numpy.nan # empty array to store data 
arousal_scale_1 = visual.RatingScale(win, labels = ["Kalm", "Opgewonden"], low=1, high=9, 
                        markerStart=5, mouseOnly = True, scale = None, acceptText = "Bevestigen", 
                        showValue = False , acceptPreText = "", pos = (0.0, 0.25), 
                        acceptSize = 1.5, stretch = 2.5)
valence_scale_1 = visual.RatingScale(win, labels = ["Verdrietig", "Vrolijk"], low=1, high=9, 
                        markerStart=5, mouseOnly = True, scale = None, acceptText = "Bevestigen", 
                        showValue = False , acceptPreText = "", pos = (0.0, -0.70), 
                        acceptSize = 1.5, stretch = 2.5)
arousal_scale_2 = visual.RatingScale(win, labels = ["Kalm", "Opgewonden"], low=1, high=9, 
                        markerStart=5, mouseOnly = True, scale = None, acceptText = "Bevestigen", 
                        showValue = False , acceptPreText = "", pos = (0.0, 0.25), 
                        acceptSize = 1.5, stretch = 2.5)
valence_scale_2 = visual.RatingScale(win, labels = ["Verdrietig", "Vrolijk"], low=1, high=9, 
                        markerStart=5, mouseOnly = True, scale = None, acceptText = "Bevestigen", 
                        showValue = False , acceptPreText = "", pos = (0.0, -0.70), 
                        acceptSize = 1.5, stretch = 2.5)
valence = visual.ImageStim(win, image = "valence.jpg", pos = (0.0, -0.45), size = (1.75,0.5))
arousal = visual.ImageStim(win, image = "arousal.jpg", pos = (0.0, 0.50), size = (1.75,0.5))
# videos randomizations 
videos_list = numpy.array(["Underworld_english_trimmed.mp4","New_York_english_trimmed.mp4"])
videos_1 = numpy.array([0,1])
videos_2 = numpy.array([1,0])

# texts
Welcome_text = "Welkom bij dit experiment \n \n Druk op SPATIE om verder te gaan"
Instructions_text = ("Tijdens dit experiment zal u 2 korte video’s te zien krijgen. " + 
                    "Uw hartslag zal gemeten worden bij het bekijken van de video’s " + 
                    "met de 2 elektrodes die de onderzoeker reeds op uw borstkast " + 
                    "heeft geplakt. Voor het begin van elke video zal er een korte " +
                    "pauze gegeven worden van 5 seconden, hierna zal de video automatisch " + 
                    "afpellen. Na elke video zal u uw emotionele response moeten " + 
                    "doorgeven via een reeks afbeeldingen zoals die hieronder afgebeeld. " + 
                    "\n Als u nog vragen heeft aarzel niet om ze te stellen aan de onderzoeker." + 
                    "\n Druk op SPATIE om het experiment te beginnen")
Goodbye_text = "Het experiment is nu afgelopen \n Bedankt voor uw participatie!"
Scale_text = "Hoe voel je je na het bekijken van de video?"


## functions 
# show a text on the screen
def message(message_text= "", response_key = "space", duration = 0, height = 0.05, pos = (0.0, 0.0), color = "white"): 
    MessageOnSCreen.text    = message_text
    MessageOnSCreen.text    = message_text
    MessageOnSCreen.height  = height
    MessageOnSCreen.pos     = pos
    MessageOnSCreen.color   = color
    
    #clean the events for keyboard  
    event.clearEvents(eventType = "keyboard")
    
    MessageOnSCreen.draw()
    win.flip()
    
    if duration == 0:
        event.waitKeys(keyList = response_key)
    else:
        time.sleep(duration)

# show video 
def show_video(file = "", size = (1050,600), units = "pix", flipVert=False, 
                flipHoriz = False, loop= False, noAudio=False, 
                volume = 0.2, autostart = True):
    
    # set the stimulus
    video.filename  = file
    video.size      = size
    video.units     = units
    video.flipVert  = flipVert
    video.flipHoriz = flipHoriz
    video.loop      = loop
    video.noAudio   = noAudio
    video.volume    = volume 
    video.autoStart = autostart
    
    while video.status != constants.FINISHED:
        # draw the movie
        video.draw()
        # flip buffers so they appear on the window
        win.flip()

## Data frame
trials[:,0] = [1,2] #save block 
trials[:,6] = participants_number # save participant number
if gender == "Male": # save gender (coded)
    trials[:,7] = 0
elif gender == "Female":
    trials[:,7] = 1
elif gender == "X":
    trials[:,7] = 2
trials[:,8] = age # save age 

# select and save randomization depending on participant number (0 is erotic, 1 is neutral)
if participants_number % 2: #if remainder is 1 aka uneven 
    trials[:,1] = numpy.transpose(videos_1)
else: # if participant number is even
    trials[:,1] = numpy.transpose(videos_2)

# creating pandas dataframe from numpy array
#trialsDF = pandas.DataFrame.from_records(trials)
# add the column names
#trialsDF.columns = ["block","condition","begin_time","end_time","valence","arousal","participant_number", "gender", "age"]


## experiment 
# welcome text 
message(message_text= Welcome_text, height = 0.1)

# instructions 
example_rating.draw()
message(message_text= Instructions_text, pos = (0.0, 0.7), height = 0.05)

trialcounter = 0
while trialcounter < trials.shape[0]: 
    # show which video
    message(message_text = "video {0}".format(int(trials[trialcounter,0])), duration = 2, height = 0.1)
    
    if speedy: 
        pass 
    else: 
        #get time stamp for begin of the video
        trials[trialcounter,2] = timer.getTime()
        # show the video
        show_video(file = videos_list[int(trials[trialcounter,1])])
        # stop the movie, this frees resources too
        video.stop()
        ##get time stamp for end of the video
        trials[trialcounter,3] = timer.getTime()
    
    # text for the rating
    MessageOnSCreen.text = Scale_text
    MessageOnSCreen.pos  = (0.0, 0.85)
    MessageOnSCreen.height  = 0.075
      
    # scales 
    if trialcounter < 1: 
        arousal_scale = arousal_scale_1
        valence_scale = valence_scale_1
    elif trialcounter == 1: 
        arousal_scale = arousal_scale_2
        valence_scale = valence_scale_2

    while arousal_scale.noResponse or valence_scale.noResponse:
        MessageOnSCreen.draw()
        valence.draw()
        arousal.draw()
        arousal_scale.draw()
        valence_scale.draw()
        
        win.flip()
    
    #save rating
    trials[trialcounter,5] = arousal_scale.getRating()
    trials[trialcounter,4] = valence_scale.getRating()
     
    trialcounter = trialcounter + 1
    
    
print(trials)

# Export as a comma separated file
numpy.savetxt(file_name, trials, fmt='%d', delimiter = ',', header= 'block, condition, valence, arousal, participant_number, gender, age') 

# end of the experiment 
message(message_text = Goodbye_text, duration = 2, height = 0.1)
win.close()
core.quit()

