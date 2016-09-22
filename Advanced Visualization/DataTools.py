# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:41:55 2014

@author: p04102916
"""

import numpy as np
import matplotlib.pyplot as plt

def dataselector(filename):
    global ismember
    global sub_set
    global sub_set_ind
    global txt

    data = np.genfromtxt(filename, delimiter=",")
    
    cc = np.shape(data)[1]
    
    datafull = np.empty((np.shape(data)[0],10))
    datafull[:,:] = np.NaN
    datafullO = np.empty((np.shape(data)[0],10))

    fig = plt.figure()
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=1)
    ax2 = plt.subplot2grid((3,3), (0,1), colspan=2)
    ax3 = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2)
    
    if cc < 10:
        for i in range(0,cc):
            datafull[:,i] = data[:,i]
    else:
        datafull = data
    
    indexes = np.arange(0, np.shape(data)[0], 1.0)
    datatmp = np.roll(datafull,1)
    datafull = datatmp
    datafull[:,0] = indexes
    
    datafullO[:,:] = datafull[:,:]

    for i in range(1,9):
        minval = min(datafull[:,i])
        if minval < 0:
            datafull[:,i] = datafull[:,i] + abs(minval)
        else:
            datafull[:,i] = datafull[:,i] - minval        
        if max(datafull[:,i]) == 0:
            datamax = 1.0
        else:
            datamax = max(datafull[:,i])
        datafull[:,i] = datafull[:,i] / datamax 
    
    ax2.axes.set_ylim(-0.1, 1.3)
    ax3.axes.set_ylim(-0.1, 1.1)
    
    L1, = ax2.plot(datafull[:,0], datafull[:,1], picker=5)
    L2, = ax2.plot(datafull[:,0], datafull[:,2], picker=5)
    L3, = ax2.plot(datafull[:,0], datafull[:,3], picker=5)
    L4, = ax2.plot(datafull[:,0], datafull[:,4], picker=5)
    L5, = ax2.plot(datafull[:,0], datafull[:,5], picker=5)
    L6, = ax2.plot(datafull[:,0], datafull[:,6], picker=5)
    L7, = ax2.plot(datafull[:,0], datafull[:,7], picker=5)
    L8, = ax2.plot(datafull[:,0], datafull[:,8], picker=5)
    L9, = ax2.plot(datafull[:,0], datafull[:,9], picker=5)
    Ls1, = ax3.plot(datafull[:,0], datafull[:,1], visible=False, picker=5)
    Ls2, = ax3.plot(datafull[:,0], datafull[:,2], visible=False, picker=5)
    Ls3, = ax3.plot(datafull[:,0], datafull[:,3], visible=False, picker=5)
    Ls4, = ax3.plot(datafull[:,0], datafull[:,4], visible=False, picker=5)
    Ls5, = ax3.plot(datafull[:,0], datafull[:,5], visible=False, picker=5)
    Ls6, = ax3.plot(datafull[:,0], datafull[:,6], visible=False, picker=5)
    Ls7, = ax3.plot(datafull[:,0], datafull[:,7], visible=False, picker=5)
    Ls8, = ax3.plot(datafull[:,0], datafull[:,8], visible=False, picker=5)
    Ls9, = ax3.plot(datafull[:,0], datafull[:,9], visible=False, picker=5)
    
    x = [0.6, 1.6]
    y = [1.5, 1.5]
    
    ax1.plot(x, y, 'bo', picker=5, ms=10)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax1.axes.set_xlim(0, 2)
    ax1.axes.set_ylim(0, 6)
    
    sub_set_ind = 0
    sub_set = 1
    
    ax1.text(0.1, 4.9, 'Current subset', fontsize=12)
    ax1.text(0.05, 2.5, 'Save & exit', fontsize=12)
    ax1.text(1.4, 2.5, 'Next', fontsize=12)
    txt = ax1.text(1.6, 4.9, sub_set, fontsize=12)
    
    ax2.text(np.round(np.shape(data)[0]/2,0), 1.1, 'Full data set', fontsize=12, ha='center')
    ax3.text(np.round(np.shape(data)[0]/2,0), 1.0, 'Data subset', fontsize=12, ha='center')
    
    #def onclick(event):
    #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #            event.button, event.x, event.y, event.xdata, event.ydata)
    #    print event.inaxes
    
    ismember = np.zeros((5,10))
#    ismember[:,0] = 1
        
    def onpick(event):
        global ismember
        global sub_set
        global sub_set_ind
        global txt
        mouseevent = event.mouseevent
        axs = mouseevent.inaxes
        thisline = event.artist
        ind = event.ind
        if axs == ax1:
            if ind == 0:
                for i in range(0,sub_set_ind+1):
                    if np.sum(ismember[i,:]) > 1:
                        selection = np.nonzero(ismember[i,:])
                        tmp1 = selection[0][:]
                        tmp2 = tuple(tmp1)
                        newdata = datafullO[:,tmp2]
                        filenamen = filename[0:filename.find('.')] + '_S' + str(i+1) + filename[filename.find('.'):len(filename)]
                        np.savetxt(filenamen, newdata, fmt='%1.5f', delimiter=",")
                plt.close()
            elif ind == 1:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                sub_set_ind = sub_set_ind + 1
                sub_set = sub_set + 1
                txt.remove()
                txt = ax1.text(1.6, 4.9, sub_set, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
        elif axs == ax2:
            # reset all lines to thin
    #        thisline.set_lw(5) # make selected line thick
    #        ax2.figure.canvas.draw() # this line is critical to change the linewidth
            if thisline == L1:
                Ls1.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,1] = 1
                thisline.set_lw(1)
                ax2.figure.canvas.draw() # this line is critical to change the linewidth
            elif thisline == L2:
                Ls2.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,2] = 1
            elif thisline == L3:
                Ls3.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,3] = 1
            elif thisline == L4:
                Ls4.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,4] = 1
            elif thisline == L5:
                Ls5.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,5] = 1
            elif thisline == L6:
                Ls6.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,6] = 1
            elif thisline == L7:
                Ls7.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,7] = 1
            elif thisline == L8:
                Ls8.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,8] = 1
            elif thisline == L9:
                Ls9.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,9] = 1
        elif axs == ax3:
            if thisline == Ls1:
                Ls1.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,1] = 0
            elif thisline == Ls2:
                Ls2.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,2] = 0
            elif thisline == Ls3:
                Ls3.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,3] = 0
            elif thisline == Ls4:
                Ls4.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,4] = 0
            elif thisline == Ls5:
                Ls5.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,5] = 0
            elif thisline == Ls6:
                Ls6.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,6] = 0
            elif thisline == Ls7:
                Ls7.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,7] = 0
            elif thisline == Ls8:
                Ls8.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,8] = 0
            elif thisline == Ls9:
                Ls9.set_visible(False)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
                ismember[sub_set_ind,9] = 0
            
    #    ax1.figure.canvas.draw() # this line is critical to change the linewidth
    
    #cid = fig.canvas.mpl_connect('button_press_event', onclick)
    cid = fig.canvas.mpl_connect('pick_event', onpick)

def dataslicer(filename):

    global ind_max
    global slices
    global txt
    global left_exist
    global leftgrey
    global left_pos
    global right_exist
    global rightgrey
    global right_pos
    global first_click
    global last_click
    global curr_slice
    global cc

    data = np.genfromtxt(filename, delimiter=",")
    
    max_slices = 5
    
    cc = np.shape(data)[1]
    
    datafull = np.empty((np.shape(data)[0],10))
    datafull[:,:] = np.NaN
    datafullO = np.empty((np.shape(data)[0],10))
    
    fig = plt.figure()
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=1)
    ax2 = plt.subplot2grid((3,3), (0,1), colspan=2)
    ax3 = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan=2)
    
    if cc < 10:
        for i in range(0,cc):
            datafull[:,i] = data[:,i]
    else:
        datafull = data
    
    datafullO[:,:] = datafull[:,:]
    
    indexes = np.arange(0, np.shape(data)[0], 1.0)
    datatmp = np.roll(datafull,1)
    datafull = datatmp
    datafull[:,0] = indexes
    
    for i in range(1,9):
        minval = min(datafull[:,i])
        if minval < 0:
            datafull[:,i] = datafull[:,i] + abs(minval)
        else:
            datafull[:,i] = datafull[:,i] - minval        
        if max(datafull[:,i]) == 0:
            datamax = 1.0
        else:
            datamax = max(datafull[:,i])
        datafull[:,i] = datafull[:,i] / datamax 
    
    ax2.axes.set_ylim(-0.1, 1.3)
    ax3.axes.set_ylim(-0.1, 1.1)
    
    L1, = ax2.plot(datafull[:,0], datafull[:,1], picker=5)
    L2, = ax2.plot(datafull[:,0], datafull[:,2], picker=5)
    L3, = ax2.plot(datafull[:,0], datafull[:,3], picker=5)
    L4, = ax2.plot(datafull[:,0], datafull[:,4], picker=5)
    L5, = ax2.plot(datafull[:,0], datafull[:,5], picker=5)
    L6, = ax2.plot(datafull[:,0], datafull[:,6], picker=5)
    L7, = ax2.plot(datafull[:,0], datafull[:,7], picker=5)
    L8, = ax2.plot(datafull[:,0], datafull[:,8], picker=5)
    L9, = ax2.plot(datafull[:,0], datafull[:,9], picker=5)
    Ls1, = ax3.plot(datafull[:,0], datafull[:,1], picker=5)
    Ls2, = ax3.plot(datafull[:,0], datafull[:,2], picker=5)
    Ls3, = ax3.plot(datafull[:,0], datafull[:,3], picker=5)
    Ls4, = ax3.plot(datafull[:,0], datafull[:,4], picker=5)
    Ls5, = ax3.plot(datafull[:,0], datafull[:,5], picker=5)
    Ls6, = ax3.plot(datafull[:,0], datafull[:,6], picker=5)
    Ls7, = ax3.plot(datafull[:,0], datafull[:,7], picker=5)
    Ls8, = ax3.plot(datafull[:,0], datafull[:,8], picker=5)
    Ls9, = ax3.plot(datafull[:,0], datafull[:,9], picker=5)
        
    x = [0.6, 1.6]
    y = [1.5, 1.5]
    
    ax1.plot(x, y, 'bo', picker=5, ms=10)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax1.axes.set_xlim(0, 2)
    ax1.axes.set_ylim(0, 6)
    
    sub_set_ind = 0
    sub_set = 1
    
    ax1.text(0.1, 4.9, 'Current slice', fontsize=12)
    ax1.text(0.05, 2.5, 'Save & exit', fontsize=12)
    ax1.text(1.4, 2.5, 'Next', fontsize=12)
    txt = ax1.text(1.6, 4.9, sub_set, fontsize=12)
    
    ax2.text(np.round(np.shape(data)[0]/2,0), 1.1, 'Full data set', fontsize=12, ha='center')
    #ax3.text(np.round(np.shape(data)[0]/2,0), 1.0, 'Data subset', fontsize=12, ha='center')
    
    #def onclick(event):
    #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #            event.button, event.x, event.y, event.xdata, event.ydata)
    #    print event.inaxes
    
    #ismember = np.zeros((5,10))
    #ismember[:,0] = 1
    slices = np.zeros((max_slices,2))
    ind_max = datafull[np.shape(datafull)[0]-1,0]
    left_exist = 0
    right_exist = 0
    left_pos = 0
    right_pos = ind_max
    first_click = 0
    last_click = 'left'
    zooming = 0
    curr_slice = 1
    
    ax2.axes.set_xlim(0, ind_max)         
    ax3.axes.set_xlim(0, ind_max)         
    
       
    def onpick(event):
        global ind_max
        global slices
        global txt
        global left_exist
        global leftgrey
        global left_pos
        global right_exist
        global rightgrey
        global right_pos
        global first_click
        global last_click
        global curr_slice
        global cc
        mouseevent = event.mouseevent
        axs = mouseevent.inaxes
    #    thisline = event.artist
        ind = event.ind
        if axs == ax1:
            if ind == 0:
                if curr_slice == 1:
    #                index_low = np.where(datafull[:,0]==round(left_pos));
    #                index_high = np.where(datafull[:,0]==round(right_pos));
                    index_low = (np.abs(datafull[:,0]-left_pos)).argmin()
                    index_high = (np.abs(datafull[:,0]-right_pos)).argmin()
                    newdata = datafullO[index_low:index_high,0:cc]
    #                newdata[:,0] = newdata[:,0] - newdata[0,0]
                    filenamen = filename[0:filename.find('.')] + '_SL1' + filename[filename.find('.'):len(filename)]
                    np.savetxt(filenamen, newdata, fmt='%1.5f', delimiter=",")
                else:
                    slices[curr_slice-1,0] = left_pos
                    slices[curr_slice-1,1] = right_pos
                    for i in range(1,curr_slice+1):
                        index_low = (np.abs(datafull[:,0]-slices[i-1,0])).argmin()
                        index_high = (np.abs(datafull[:,0]-slices[i-1,1])).argmin()
                        newdata = datafullO[index_low:index_high,0:cc]
    #                    newdata[:,0] = newdata[:,0] - newdata[0,0]
                        filenamen = filename[0:filename.find('.')] + '_SL' + str(i) + filename[filename.find('.'):len(filename)]
                        np.savetxt(filenamen, newdata, fmt='%1.5f', delimiter=",")
    #                if np.sum(ismember[i,:]) > 1:
    #                    selection = np.nonzero(ismember[i,:])
    #                    tmp1 = selection[0][:]
    #                    tmp2 = tuple(tmp1)
    #                    newdata = data[:,tmp2]
                plt.close()
            elif ind == 1:
    #            for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
    #                line.set_visible(False)
    #            sub_set_ind = sub_set_ind + 1
    #            sub_set = sub_set + 1
                slices[curr_slice-1,0] = left_pos
                slices[curr_slice-1,1] = right_pos
                curr_slice = curr_slice + 1
                txt.remove()
                txt = ax1.text(1.6, 4.9, curr_slice, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
                leftgrey.remove()
                rightgrey.remove()
                left_pos = 0
                right_pos = ind_max
                left_exist = 0
                right_exist = 0
                first_click = 0
                last_click = 'left'
    #            ax2.axes.set_ylim(-0.1, 1.1)
    #            ax2.axes.set_xlim(0, ind_max)         
                ax2.figure.canvas.draw()
                ax3.axes.set_ylim(-0.1, 1.1)
                ax3.axes.set_xlim(0, ind_max)         
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
    def onclick(event):
        global left_exist
        global leftgrey
        global left_pos
        global right_exist
        global rightgrey
        global right_pos
        global first_click
        global last_click
        global ind_max
        global zooming
    #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #        event.button, event.x, event.y, event.xdata, event.ydata)
        xpos = event.xdata
    #    ypos = event.ydata
        axs = event.inaxes
        if axs == ax2:
    #        plt.axvline(x=1) 
    #        xl = [xpos, xpos]
    #        yl = [0, 1000]
    #        ax2.axvline(xpos)
            if first_click == 0:
                leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                first_click = 1
                left_pos = xpos
                left_exist = 1
                ax2.figure.canvas.draw() # this line is critical to change the linewidth
    
                minval = 999999.0
                maxval = -999999.0
                for i in range(1,9):
                    if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                        if min(datafull[int(left_pos):int(right_pos),i])<minval:
                            minval = min(datafull[int(left_pos):int(right_pos),i])
                        if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                            maxval = max(datafull[int(left_pos):int(right_pos),i])
                ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            elif first_click == 1:
                if xpos < left_pos:
                    xpos = left_pos + 10
                rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                first_click = 2
                right_pos = xpos
                right_exist = 1
                ax2.figure.canvas.draw() # this line is critical to change the linewidth
    
                minval = 999999.0
                maxval = -999999.0
                for i in range(1,9):
                    if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                        if min(datafull[int(left_pos):int(right_pos),i])<minval:
                            minval = min(datafull[int(left_pos):int(right_pos),i])
                        if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                            maxval = max(datafull[int(left_pos):int(right_pos),i])
                ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            elif first_click == 2:
                middle = (left_pos + right_pos) / 2
                if xpos < left_pos:
                     leftgrey.remove()
                     leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                     left_pos = xpos
                     ax2.figure.canvas.draw() # this line is critical to change the linewidth
    
                     minval = 999999.0
                     maxval = -999999.0
                     for i in range(1,9):
                         if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                             if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                 minval = min(datafull[int(left_pos):int(right_pos),i])
                             if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                 maxval = max(datafull[int(left_pos):int(right_pos),i])
                     ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                     ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
                elif xpos > right_pos:
                     rightgrey.remove()
                     rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                     ax2.figure.canvas.draw() # this line is critical to change the linewidth
                     right_pos = xpos
    
                     minval = 999999.0
                     maxval = -999999.0
                     for i in range(1,9):
                         if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                             if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                 minval = min(datafull[int(left_pos):int(right_pos),i])
                             if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                 maxval = max(datafull[int(left_pos):int(right_pos),i])
                     ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                     ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
                elif xpos > middle:
                     rightgrey.remove()
                     rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                     ax2.figure.canvas.draw() # this line is critical to change the linewidth
                     right_pos = xpos
    
                     minval = 999999.0
                     maxval = -999999.0
                     for i in range(1,9):
                         if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                             if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                 minval = min(datafull[int(left_pos):int(right_pos),i])
                             if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                 maxval = max(datafull[int(left_pos):int(right_pos),i])
                     ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                     ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
                elif xpos < middle: 
                     leftgrey.remove()
                     leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                     left_pos = xpos
                     ax2.figure.canvas.draw() # this line is critical to change the linewidth
    
                     minval = 999999.0
                     maxval = -999999.0
                     for i in range(1,9):
                         if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                             if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                 minval = min(datafull[int(left_pos):int(right_pos),i])
                             if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                 maxval = max(datafull[int(left_pos):int(right_pos),i])
                     ax3.axes.set_ylim(minval-0.1, maxval+0.1)
                     ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            ax3.axes.set_xlim(left_pos, right_pos)         
    #        ax3.axes.set_ylim(-0.1, 1.1)
            ax3.figure.canvas.draw() # this line is critical to change the linewidth
        elif axs == ax3:
            if first_click == 2:
                middle = (left_pos + right_pos) / 2
                if xpos < middle:
                    leftgrey.remove()
                    leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                    left_pos = xpos
                    ax2.figure.canvas.draw() # this line is critical to change the linewidth

                    minval = 999999.0
                    maxval = -999999.0
                    for i in range(1,9):
                        if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                            if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                minval = min(datafull[int(left_pos):int(right_pos),i])
                            if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                maxval = max(datafull[int(left_pos):int(right_pos),i])
                    ax3.axes.set_ylim(minval-0.1, maxval+0.1)

                elif xpos > middle:
                    rightgrey.remove()
                    rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                    right_pos = xpos
                    ax2.figure.canvas.draw() # this line is critical to change the linewidth

                    minval = 999999.0
                    maxval = -999999.0
                    for i in range(1,9):
                        if str(min(datafull[int(left_pos):int(right_pos),i]))!='nan':
                            if min(datafull[int(left_pos):int(right_pos),i])<minval:
                                minval = min(datafull[int(left_pos):int(right_pos),i])
                            if max(datafull[int(left_pos):int(right_pos),i])>maxval:
                                maxval = max(datafull[int(left_pos):int(right_pos),i])
                    ax3.axes.set_ylim(minval-0.1, maxval+0.1)

                ax3.axes.set_xlim(left_pos, right_pos)         
#                ax3.axes.set_ylim(-0.1, 1.1)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            
    cid = fig.canvas.mpl_connect('pick_event', onpick)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    




