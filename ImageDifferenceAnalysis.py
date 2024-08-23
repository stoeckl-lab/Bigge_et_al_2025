# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 16:17:43 2020

@author: rob38dk
"""

#%% Select Video
import tkinter as tk
from tkinter import filedialog

root=tk.Tk()
root.attributes("-topmost", True)
root.withdraw()
file_path=filedialog.askopenfilename(initialdir="X:\ImageDifference\Videos", title="Select a Video")
print(file_path)
condition= input('Enter Condition name: ')


if "turn" in file_path:
    file_name=file_path[-41:]
else:    
    file_name= file_path[-36:]

#%% Downsample Video
import numpy as np
import cv2 as cv
import os.path

count1=0
count2=0

cap = cv.VideoCapture(file_path)
frame_width = int(cap.get(3))
frame_height = int(cap.get(4))

new_dir= os.path.join("X:\ImageDifference\Analysis", condition[:-2], condition)
if not os.path.exists(new_dir):
    os.makedirs(new_dir)
    
os.chdir(new_dir)
out = cv.VideoWriter('downsample_'+file_name,cv.VideoWriter_fourcc('M','J','P','G'), 10, (frame_width,frame_height))

while(1):
    ret, frame = cap.read()
    if ret == False:
        break
    
    if count1 == count2:
        out.write(frame)
        count2=count2 +20
        cv.imshow('frame',frame)
        if cv.waitKey(1) & 0xFF == ord('q'):
            break
        
    count1= count1+1    
    
cap.release()
out.release()


#%% Blur image with gaussian filter and save as new video
#import scipy.ndimage as ndimage

cap=cv.VideoCapture('downsample_'+file_name)
blur= cv.VideoWriter('down_blur_'+file_name,cv.VideoWriter_fourcc('M','J','P','G'), 10, (frame_width,frame_height))


while(1):
    ret,fr = cap.read()
    if ret == False:
        break
    
    #gaussimg= ndimage.gaussian_filter(fr, sigma=(5,5,0), order=0)
    #img=cv.imread(fr)
    
    gaussimg= cv.GaussianBlur(fr,(5,5),0)
    #cv.imshow(gaussimg)
    blur.write(gaussimg)
    #cv.imshow('Blurr',gaussimg)
    #if cv.waitkey(1) & 0xFF == ord('q'):
        #break

cap.release()
blur.release()    


#%% Calculate Optical Flow of the downsampled Video

cap = cv.VideoCapture(cv.samples.findFile('down_blur_'+file_name))
ret, frame1 = cap.read()
prvs = cv.cvtColor(frame1,cv.COLOR_BGR2GRAY)
hsv = np.zeros_like(frame1)
hsv[...,1] = 255
angle=[]
magnitude=[]

while(1):
    ret, frame2 = cap.read()
    if ret == False:
        break
    nex = cv.cvtColor(frame2,cv.COLOR_BGR2GRAY)
    flow = cv.calcOpticalFlowFarneback(prvs,nex, None, 0.5, 3, 15, 3, 5, 1.2, 0)
    mag, ang = cv.cartToPolar(flow[...,0], flow[...,1])
    hsv[...,0] = ang*180/np.pi/2
    hsv[...,2] = cv.normalize(mag,None,0,255,cv.NORM_MINMAX)
    bgr = cv.cvtColor(hsv,cv.COLOR_HSV2BGR)
    cv.imshow('frame2',bgr)
    k = cv.waitKey(30) & 0xff
    if k == 27:
        break
    # elif k == ord('s'):
    #     cv.imwrite('opticalfb.png',frame2)
    #     cv.imwrite('opticalhsv.png',bgr)
    prvs = nex
    angle.append(ang)
    magnitude.append(mag) 
    
#%% Save Angle and Magnitude for Matlab ( ??Decide if mean or median??)

# median optical Flow
m_angle= np.median(angle, axis=0)
m_magnitude=np.median(magnitude, axis=0)
    
# mean optical flow
#m_angle= np.mean(angle, axis=0)
#m_magnitude=np.mean(magnitude, axis=0)
    
    
import scipy.io
scipy.io.savemat('angle_'+file_name[:-9:]+'.mat', mdict={'angle': m_angle})
scipy.io.savemat('magnitude_'+file_name[:-9:]+'.mat', mdict={'magnitude': m_magnitude})   

AngleStack= np.empty((len(angle),),dtype=np.object)
for i in range(len(angle)):
    AngleStack[i]=angle[i]
    
scipy.io.savemat('AngleStack_'+file_name[:-9]+'.mat', {'AngleStack':AngleStack}) 

MagniStack= np.empty((len(magnitude),),dtype=np.object)
for i in range(len(magnitude)):
    MagniStack[i]=magnitude[i]
    
scipy.io.savemat('MagniStack_'+file_name[:-9]+'.mat',{'MagniStack':MagniStack})    


    