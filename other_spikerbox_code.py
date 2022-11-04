import serial
import numpy as np
import matplotlib.pyplot as plt
import time
#matplotlib notebook

def read_arduino(ser,inputBufferSize):
#    data = ser.readline(inputBufferSize)
    data = ser.read(inputBufferSize)
    out =[(int(data[i])) for i in range(0,len(data))]
    return out

def process_data(data):
    data_in = np.array(data)
    result = []
    i = 1
    while i < len(data_in)-1:
        if data_in[i] > 127:
            # Found beginning of frame
            # Extract one sample from 2 bytes
            intout = (np.bitwise_and(data_in[i],127))*128
            i = i + 1
            intout = intout + data_in[i]
            result = np.append(result,intout)
        i=i+1
    return result

baudrate = 230400
cport = 'COM6'  # set the correct port before you run it
ser = serial.Serial(port=cport, baudrate=baudrate)    

# take continuous data stream 
inputBufferSize = 10000 # keep between 2000-20000
ser.timeout = inputBufferSize/20000.0  # set read timeout, 20000 is one second
ser.set_buffer_size(rx_size = inputBufferSize)

total_time = 5.0; # time in seconds [[1 s = 20000 buffer size]] // 8 MINUTES (480 s) FOR TESTING 
max_time = 5.0; # time plotted in window [s]
N_loops = 20000.0/inputBufferSize*total_time

T_acquire = inputBufferSize/20000.0    # length of time that data is acquired for 
N_max_loops = max_time/T_acquire    # total number of loops to cover desire time window

for k in range(0,int(N_loops)):
    data = read_arduino(ser,inputBufferSize)
    data_temp = process_data(data)
    if k <= N_max_loops:
        if k==0:
            data_plot = data_temp
        else:
            data_plot = np.append(data_plot, data_temp) #original goes (data_temp, data_plot) but data_temp is what should be appended to te already existing data_plot
        t = (min(k+1,N_max_loops))*inputBufferSize/20000.0*np.linspace(0,1,(data_plot).size)
    else:
        data_plot = np.roll(data_plot,len(data_temp))
        data_plot[0:len(data_temp)] = data_temp
    t = (min(k+1,N_max_loops))*inputBufferSize/20000.0*np.linspace(0,1,(data_plot).size)

    
# save the plot above when ready
name_of_file = 'example_file.txt'
np.savetxt(name_of_file, np.c_[t, np.real(data_plot)])
