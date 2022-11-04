# Backyard Brains Sep. 2019
# Made for python 3
# First install serial library
# Install numpy, pyserial, matplotlib
# pip3 install pyserial
#
# Code will read, parse and display data from BackyardBrains' serial devices
#
# Written by Stanislav Mircic
# stanislav@backyardbrains.com

import threading
import serial
import time
import matplotlib.pyplot as plt 
import numpy as np
#global makes a variabel global and thus it exist outside of a function. Not really necessary since this is alreayd outside a funciton
global connected 
connected = False
#change name of the port here
#port = 'COM4' #check in device manager (windows name for port) 
port = '/dev/cu.usbmodem1411101'# linux version of port name 
baud = 230400 #How fast your COM port operates. this is a standard value (among other options)
global input_buffer
global sample_buffer
global cBufTail
cBufTail = 0 #variable is equal to 0
input_buffer = [] #variable is an empty list 
sample_rate = 10000
display_size = 30000 #3 seconds (can be changed)
sample_buffer = np.linspace(0,0,display_size) #numpy array from 0 to 0 with display_size elements 
serial_port = serial.Serial(port, baud, timeout=0) # read the port with a specific baud and return immediately in any case, returning zero or more, up to the requested number of bytes




def checkIfNextByteExist():
        global cBufTail
        global input_buffer
        tempTail = cBufTail + 1
        
        if tempTail==len(input_buffer): 
            return False
        return True
    

def checkIfHaveWholeFrame():
        global cBufTail
        global input_buffer
        tempTail = cBufTail + 1
        while tempTail!=len(input_buffer): 
            nextByte  = input_buffer[tempTail] & 0xFF #mask the first bytes and leave just the last one or last 8 (not clear)
            if nextByte > 127:
                return True
            tempTail = tempTail +1
        return False;
    
def areWeAtTheEndOfFrame():
        global cBufTail
        global input_buffer
        tempTail = cBufTail + 1
        nextByte  = input_buffer[tempTail] & 0xFF
        if nextByte > 127:
            return True
        return False

def numberOfChannels():
    return 1

def handle_data(data):
    global input_buffer
    global cBufTail
    global sample_buffer    
    if len(data)>0:

        cBufTail = 0
        haveData = True
        weAlreadyProcessedBeginingOfTheFrame = False
        numberOfParsedChannels = 0
        
        while haveData:
            MSB  = input_buffer[cBufTail] & 0xFF
            
            if(MSB > 127):
                weAlreadyProcessedBeginingOfTheFrame = False
                numberOfParsedChannels = 0 #increases 1 with each reiteration 
                
                if checkIfHaveWholeFrame():
                    
                    while True: #the loop will be infinate, hence the break if loops
                        
                        MSB  = input_buffer[cBufTail] & 0xFF
                        if(weAlreadyProcessedBeginingOfTheFrame and (MSB>127)): #break from the while loop
                            #we have begining of the frame inside frame
                            #something is wrong
                            break #continue as if we have new frame
            
                        MSB  = input_buffer[cBufTail] & 0x7F
                        weAlreadyProcessedBeginingOfTheFrame = True
                        cBufTail = cBufTail +1
                        LSB  = input_buffer[cBufTail] & 0xFF

                        if LSB>127: #break from the while loop
                            break #continue as if we have new frame

                        LSB  = input_buffer[cBufTail] & 0x7F
                        MSB = MSB<<7
                        writeInteger = LSB | MSB #set to 1 is one of the bytes is 1
                        numberOfParsedChannels = numberOfParsedChannels+1
                        if numberOfParsedChannels>numberOfChannels(): #break from the while loop. numberOfChannels is 1 
            
                            #we have more data in frame than we need
                            #something is wrong with this frame
                            break #continue as if we have new frame
            

                        sample_buffer = np.append(sample_buffer,writeInteger-512) #why 512
                        

                        if areWeAtTheEndOfFrame(): #break from the while loop, were at the end of the frame
                            break
                        else:
                            cBufTail = cBufTail +1
                else:
                    haveData = False
                    break
            if(not haveData): #become if true if haveData is false. end up in the while loop with haveData
                break
            cBufTail = cBufTail +1
            if cBufTail==len(input_buffer):
                haveData = False
                break


def read_from_port(ser): #Ser becomes the serial_port. I assume this is recording (BUT NOT SAVING) the data from the spikerbox 
    global connected
    global input_buffer
    while not connected:
        #serin = ser.read()
        connected = True

        while True:
           
           reading = ser.read(1024) #read 1024 bytes at a time (?)
           if(len(reading)>0):
                reading = list(reading)
#here we overwrite if we left some parts of the frame from previous processing 
#should be changed             
                input_buffer = reading.copy()
                print("len(reading)",len(reading))
                handle_data(reading) #all functions are created to be able to write this line 
           
           time.sleep(0.001)

thread = threading.Thread(target=read_from_port, args=(serial_port,)) #target is invoked by run(), args is the argument needed for the target 
thread.start()
xi = np.linspace(-display_size/sample_rate, 0, num=display_size) # from - number of seconds to 0, get display size amount of samples which are equally spaced 
#xi is the time variable 
while True: #infite loop 
    plt.ion() #interactive mode is on 
    plt.show(block=False) #all figure windows are displayed and return immediately
    if(len(sample_buffer)>0):
        #i = len(sample_buffer)
        print(len(sample_buffer))
        yi = sample_buffer.copy()
        yi = yi[-display_size:]
        sample_buffer = sample_buffer[-display_size:]
        plt.clf() # clear the current figure 

        plt.ylim(-550, 550)
        plt.plot(xi, yi, linewidth=1, color='royalblue')
        plt.pause(0.001)
        time.sleep(0.08)
        