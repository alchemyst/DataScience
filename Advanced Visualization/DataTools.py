# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 10:41:55 2014

@author: p04102916
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_lines(ax, datafull, N=10, **kwargs):
    return [ax.plot(datafull[:, 0], datafull[:, i], picker=5, **kwargs)[0]
            for i in range(1, N)]

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

    L1, L2, L3, L4, L5, L6, L7, L8, L9 = plot_lines(ax2, datafull)
    Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9 = plot_lines(ax3, datafull, visible=False)

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

        if (fig.canvas.manager.toolbar._active is None) == True:

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
    
    L1, L2, L3, L4, L5, L6, L7, L8, L9 = plot_lines(ax2, datafull)
    Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9 = plot_lines(ax3, datafull)
    
    Ls0, = ax3.plot([np.round(np.shape(data)[0]/2,0), np.round(np.shape(data)[0]/2,0)], [-0.1, 1.5], 'k', linestyle='dashed')
       
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

        if (fig.canvas.manager.toolbar._active is None) == True:

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
        
                    Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
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
        
                    Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
                    ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                if first_click == 2:
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
        
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
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
        
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
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
        
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
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
        
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
        
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
    
                        Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
    
                        ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
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
    
                        Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
    
                        ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
                    ax3.axes.set_xlim(left_pos, right_pos)         
    #                ax3.axes.set_ylim(-0.1, 1.1)
                    ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            
    cid = fig.canvas.mpl_connect('pick_event', onpick)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    
def datamunge(filename):

    global curr_select
    global curr_task
    global x_pos1
    global x_pos2
    global filtO
    global txtf

    data = np.genfromtxt(filename, delimiter=",")
    
    cc = np.shape(data)[1]
    cr = np.shape(data)[0]
    
    datafull = np.empty((cr,10))
    datafull[:,:] = np.NaN
    datafullO = np.empty((cr,10))
    datafullScO = np.empty((cr,10))
    dataOut = np.empty((cr,cc))
    
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
    
    dataSc = np.empty((9,2))
    
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
        dataSc[i-1,0] = minval
        dataSc[i-1,1] = datamax
    
    datafullScO[:,:] = datafull[:,:]
    
    
    ax2.axes.set_ylim(-0.1, 1.3)
    ax3.axes.set_ylim(-0.1, 1.1)
    
    L1, L2, L3, L4, L5, L6, L7, L8, L9 = plot_lines(ax2, datafull)
    Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9 = plot_lines(ax3, datafull, visible=False)
    L1o, L2o, L3o, L4o, L5o, L6o, L7o, L8o, L9o = plot_lines(ax3, datafull, lw=4, alpha=0.2)
    Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o = plot_lines(ax3, datafull, lw=4, alpha=0.2, visible=False)
    
    x = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8]
    y = [0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 0.6, 1.6, 2.6, 3.6, 4.6, 5.6]
    #x = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.8, 1.8, 1.8, 1.8]
    #y = [0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 0.6, 1.6, 4.6, 5.6]
    
    ax1.plot(x, y, 'bo', picker=5, ms=10)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax1.axes.set_xlim(0, 2)
    ax1.axes.set_ylim(0, 6.1)
    
    sub_set_ind = 0
    sub_set = 1
        
    Ls0, = ax1.plot([1.0, 1.0], [-0.1, 6.5], 'k')
    
    ax1.text(0.4, 0.35, 'HLG L', fontsize=12)
    ax1.text(0.4, 1.35, 'HLG R', fontsize=12)
    ax1.text(0.4, 2.35, 'Interp', fontsize=12)
    ax1.text(0.4, 3.35, 'CAV', fontsize=12)
    ax1.text(0.4, 4.35, 'RlinD', fontsize=12)
    ax1.text(1.1, 1.35, 'Reset', fontsize=12)
    ax1.text(1.1, 0.35, 'Exit', fontsize=12)
    #txt = ax1.text(1.6, 4.9, sub_set, fontsize=12)
    
    filtBx = [1.8, 1.8]
    filtBy = [3.6, 2.6]
    filtO = 1
    
    txtf = ax1.text(1.1, 2.8, filtO, fontsize=15)
    ax1.plot(filtBx[0], filtBy[0], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx[1], filtBy[1], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx[0], filtBy[0], 'b^', ms=8)
    ax1.plot(filtBx[1], filtBy[1], 'bv', ms=8)
    
    ax2.text(np.round(np.shape(data)[0]/2,0), 1.1, 'Full data set', fontsize=12, ha='center')
    ax3.text(np.round(np.shape(data)[0]/2,0), 1.0, 'Selected input', fontsize=12, ha='center')
    
    #def onclick(event):
    #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
    #            event.button, event.x, event.y, event.xdata, event.ydata)
    #    print event.inaxes
    
    ismember = np.zeros((5,10))
    #ismember[:,0] = 1
    
    curr_task = 0
    curr_select = 0
    x_pos1 = -9999.99
    x_pos2 = -9999.99
        
    def onpick(event):
    #    global ismember
    #    global sub_set
    #    global sub_set_ind
    #    global txt
    
        global curr_select
        global curr_task
        global x_pos1
        global x_pos2
        global filtO
        global txtf
    
        mouseevent = event.mouseevent
        axs = mouseevent.inaxes
        thisline = event.artist
        ind = event.ind
        if axs == ax1:
    #        print(ind)
            if ind == 0:
                curr_task = 1
                x_pos1 = -9999.99
                x_pos2 = -9999.99
            elif ind == 1:
                curr_task = 2
                x_pos1 = -9999.99
                x_pos2 = -9999.99
            elif ind == 2:
                curr_task = 3
                x_pos1 = -9999.99
                x_pos2 = -9999.99
            elif ind == 3:
                if (curr_select != 0):
                    def movingaverage(interval, window_size):
                        window= np.ones(int(window_size))/float(window_size)
                        return np.convolve(interval, window, 'same')
                    if Ls1.get_visible():
                        datafull[:,1] = movingaverage(datafull[:,1], filtO)
                        Ls1.set_ydata(datafull[:,1])
                        L1.set_ydata(datafull[:,1])
                    elif Ls2.get_visible():
                        datafull[:,2] = movingaverage(datafull[:,2], filtO)
                        Ls2.set_ydata(datafull[:,2])
                        L2.set_ydata(datafull[:,2])
                    elif Ls3.get_visible():
                        datafull[:,3] = movingaverage(datafull[:,3], filtO)
                        Ls3.set_ydata(datafull[:,3])
                        L3.set_ydata(datafull[:,3])
                    elif Ls4.get_visible():
                        datafull[:,4] = movingaverage(datafull[:,4], filtO)
                        Ls4.set_ydata(datafull[:,4])
                        L4.set_ydata(datafull[:,4])
                    elif Ls5.get_visible():
                        datafull[:,5] = movingaverage(datafull[:,5], filtO)
                        Ls5.set_ydata(datafull[:,5])
                        L5.set_ydata(datafull[:,5])
                    elif Ls6.get_visible():
                        datafull[:,6] = movingaverage(datafull[:,6], filtO)
                        Ls6.set_ydata(datafull[:,6])
                        L6.set_ydata(datafull[:,6])
                    elif Ls7.get_visible():
                        datafull[:,7] = movingaverage(datafull[:,7], filtO)
                        Ls7.set_ydata(datafull[:,7])
                        L7.set_ydata(datafull[:,7])
                    elif Ls8.get_visible():
                        datafull[:,8] = movingaverage(datafull[:,8], filtO)
                        Ls8.set_ydata(datafull[:,8])
                        L8.set_ydata(datafull[:,8])
                    elif Ls9.get_visible():
                        datafull[:,9] = movingaverage(datafull[:,9], filtO)
                        Ls9.set_ydata(datafull[:,9])
                        L9.set_ydata(datafull[:,9])
                    
                    ax3.figure.canvas.draw() # this line is critical to change the linewidth
                    curr_task = 0
                    x_pos1 = -9999.99
                    x_pos2 = -9999.99
    
            elif ind == 4:
                curr_task = 5
                x_pos1 = -9999.99
                x_pos2 = -9999.99
    
            elif ind == 7:
                if Ls1.get_visible():
                    datafull[:,1] = datafullScO[:,1]
                    Ls1.set_ydata(datafull[:,1])
                    L1.set_ydata(datafull[:,1])
                elif Ls2.get_visible():
                    datafull[:,2] = datafullScO[:,2]
                    Ls2.set_ydata(datafull[:,2])
                    L2.set_ydata(datafull[:,2])
                elif Ls3.get_visible():
                    datafull[:,3] = datafullScO[:,3]
                    Ls3.set_ydata(datafull[:,3])
                    L3.set_ydata(datafull[:,3])
                elif Ls4.get_visible():
                    datafull[:,4] = datafullScO[:,4]
                    Ls4.set_ydata(datafull[:,4])
                    L4.set_ydata(datafull[:,4])
                elif Ls5.get_visible():
                    datafull[:,5] = datafullScO[:,5]
                    Ls5.set_ydata(datafull[:,5])
                    L5.set_ydata(datafull[:,5])
                elif Ls6.get_visible():
                    datafull[:,6] = datafullScO[:,6]
                    Ls6.set_ydata(datafull[:,6])
                    L6.set_ydata(datafull[:,6])
                elif Ls7.get_visible():
                    datafull[:,7] = datafullScO[:,7]
                    Ls7.set_ydata(datafull[:,7])
                    L7.set_ydata(datafull[:,7])
                elif Ls8.get_visible():
                    datafull[:,8] = datafullScO[:,8]
                    Ls8.set_ydata(datafull[:,8])
                    L8.set_ydata(datafull[:,8])
                elif Ls9.get_visible():
                    datafull[:,9] = datafullScO[:,9]
                    Ls9.set_ydata(datafull[:,9])
                    L9.set_ydata(datafull[:,9])
               
                ax2.figure.canvas.draw() # this line is critical to change the linewidth
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            elif ind == 8:
                filtO = filtO - 1
                if filtO == 0:
                    filtO = 1
                txtf.remove()
                txtf = ax1.text(1.1, 2.8, filtO, fontsize=15)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
            elif ind == 9:
                filtO = filtO + 1
                txtf.remove()
                txtf = ax1.text(1.1, 2.8, filtO, fontsize=15)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
    
            elif ind == 6:
                filenamen = filename[0:filename.find('.')] + '_F' + filename[filename.find('.'):len(filename)]
                for i in range(1,cc+1):
                    dataOut[:,i-1] = (datafull[:,i] * dataSc[i-1,1]) + dataSc[i-1,0]
                np.savetxt(filenamen, dataOut, fmt='%1.5f', delimiter=",")
                plt.close()
    
        elif axs == ax2:
            # reset all lines to thin
    #        thisline.set_lw(5) # make selected line thick
    #        ax2.figure.canvas.draw() # this line is critical to change the linewidth
            if thisline == L1:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls1.set_visible(True)
                Ls1o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,1] = 1
    #            thisline.set_lw(1)
    #            ax2.figure.canvas.draw() # this line is critical to change the linewidth
                curr_select = 1
            elif thisline == L2:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls2.set_visible(True)
                Ls2o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,2] = 1
                curr_select = 2
            elif thisline == L3:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls3.set_visible(True)
                Ls3o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,3] = 1
                curr_select = 3
            elif thisline == L4:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls4.set_visible(True)
                Ls4o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,4] = 1
                curr_select = 4
            elif thisline == L5:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls5.set_visible(True)
                Ls5o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,5] = 1
                curr_select = 5
            elif thisline == L6:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls6.set_visible(True)
                Ls6o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,6] = 1
                curr_select = 6
            elif thisline == L7:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls7.set_visible(True)
                Ls7o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,7] = 1
                curr_select = 7
            elif thisline == L8:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls8.set_visible(True)
                Ls8o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,8] = 1
                curr_select = 8
            elif thisline == L9:
                for line in [Ls1, Ls2, Ls3, Ls4, Ls5, Ls6, Ls7, Ls8, Ls9]:
                    line.set_visible(False)
                for line in [Ls1o, Ls2o, Ls3o, Ls4o, Ls5o, Ls6o, Ls7o, Ls8o, Ls9o]:
                    line.set_visible(False)
                Ls9.set_visible(True)
                Ls9o.set_visible(True)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,9] = 1
                curr_select = 9
    #    elif axs == ax3:
    #        if thisline == Ls1:
    #            Ls1.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,1] = 0
    #        elif thisline == Ls2:
    #            Ls2.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,2] = 0
    #        elif thisline == Ls3:
    #            Ls3.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,3] = 0
    #        elif thisline == Ls4:
    #            Ls4.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,4] = 0
    #        elif thisline == Ls5:
    #            Ls5.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,5] = 0
    #        elif thisline == Ls6:
    #            Ls6.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,6] = 0
    #        elif thisline == Ls7:
    #            Ls7.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,7] = 0
    #        elif thisline == Ls8:
    #            Ls8.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,8] = 0
    #        elif thisline == Ls9:
    #            Ls9.set_visible(False)
    #            ax3.figure.canvas.draw() # this line is critical to change the linewidth
    #            ismember[sub_set_ind,9] = 0
    #        
    #    ax1.figure.canvas.draw() # this line is critical to change the linewidth
    
    def onclick(event):
    
        global curr_select
        global curr_task
        global x_pos1
        global x_pos2
        global filtO
        global txtf
    
        if (fig.canvas.manager.toolbar._active is None) == True:
    
        #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        #        event.button, event.x, event.y, event.xdata, event.ydata)
            xpos = event.xdata
            ypos = event.ydata
            axs = event.inaxes
    
            if axs == ax3:
                if (curr_task != 0)and(curr_select != 0):
    
                    if x_pos1 == -9999.99:
                        x_pos1 = xpos
                    elif x_pos2 == -9999.99:
                        x_pos2 = xpos
    
                        if curr_task == 3:
                            if x_pos1 > x_pos2:
                                right_ind = int(x_pos1)
                                left_ind = int(x_pos2)
                            else:    
                                right_ind = int(x_pos2)
                                left_ind = int(x_pos1)
                            if Ls1.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,1])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),1])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,1] = interp
                                Ls1.set_ydata(datafull[:,1])
                                L1.set_ydata(datafull[:,1])
                            elif Ls2.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,2])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),2])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,2] = interp
                                Ls2.set_ydata(datafull[:,2])
                                L2.set_ydata(datafull[:,2])
                            elif Ls3.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,3])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),3])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,3] = interp
                                Ls3.set_ydata(datafull[:,3])
                                L3.set_ydata(datafull[:,3])
                            elif Ls4.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,4])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),4])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,4] = interp
                                Ls4.set_ydata(datafull[:,4])
                                L4.set_ydata(datafull[:,4])
                            elif Ls5.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,5])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),5])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,5] = interp
                                Ls5.set_ydata(datafull[:,5])
                                L5.set_ydata(datafull[:,5])
                            elif Ls6.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,6])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),6])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,6] = interp
                                Ls6.set_ydata(datafull[:,6])
                                L6.set_ydata(datafull[:,6])
                            elif Ls7.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,7])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),7])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,7] = interp
                                Ls7.set_ydata(datafull[:,7])
                                L7.set_ydata(datafull[:,7])
                            elif Ls8.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,8])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),8])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,8] = interp
                                Ls8.set_ydata(datafull[:,8])
                                L8.set_ydata(datafull[:,8])
                            elif Ls9.get_visible():
                                left_val = np.mean(datafull[(left_ind-int(filtO)):left_ind,9])
                                right_val = np.mean(datafull[right_ind:(right_ind+int(filtO)),9])
                                interp = np.linspace(left_val,right_val,(right_ind-left_ind))
                                datafull[left_ind:right_ind,9] = interp
                                Ls9.set_ydata(datafull[:,9])
                                L9.set_ydata(datafull[:,9])
                            curr_task = 0
    
                        if curr_task == 5:
                            if x_pos1 > x_pos2:
                                right_ind = int(x_pos1)
                                left_ind = int(x_pos2)
                            else:    
                                right_ind = int(x_pos2)
                                left_ind = int(x_pos1)
                            if Ls1.get_visible():
                                interp = np.linspace(datafull[left_ind,1],datafull[right_ind,1],(right_ind-left_ind))
                                datafull[left_ind:right_ind,1] = datafull[left_ind:right_ind,1] - interp[:] + datafull[left_ind,1]
                                Ls1.set_ydata(datafull[:,1])
                                L1.set_ydata(datafull[:,1])
                            elif Ls2.get_visible():
                                interp = np.linspace(datafull[left_ind,2],datafull[right_ind,2],(right_ind-left_ind))
                                datafull[left_ind:right_ind,2] = datafull[left_ind:right_ind,2] - interp[:] + datafull[left_ind,2]
                                Ls2.set_ydata(datafull[:,2])
                                L2.set_ydata(datafull[:,2])
                            elif Ls3.get_visible():
                                interp = np.linspace(datafull[left_ind,3],datafull[right_ind,3],(right_ind-left_ind))
                                datafull[left_ind:right_ind,3] = datafull[left_ind:right_ind,3] - interp[:] + datafull[left_ind,3]
                                Ls3.set_ydata(datafull[:,3])
                                L3.set_ydata(datafull[:,3])
                            elif Ls4.get_visible():
                                interp = np.linspace(datafull[left_ind,4],datafull[right_ind,4],(right_ind-left_ind))
                                datafull[left_ind:right_ind,4] = datafull[left_ind:right_ind,4] - interp[:] + datafull[left_ind,4]
                                Ls4.set_ydata(datafull[:,4])
                                L4.set_ydata(datafull[:,4])
                            elif Ls5.get_visible():
                                interp = np.linspace(datafull[left_ind,5],datafull[right_ind,5],(right_ind-left_ind))
                                datafull[left_ind:right_ind,5] = datafull[left_ind:right_ind,5] - interp[:] + datafull[left_ind,5]
                                Ls5.set_ydata(datafull[:,5])
                                L5.set_ydata(datafull[:,5])
                            elif Ls6.get_visible():
                                interp = np.linspace(datafull[left_ind,6],datafull[right_ind,6],(right_ind-left_ind))
                                datafull[left_ind:right_ind,6] = datafull[left_ind:right_ind,6] - interp[:] + datafull[left_ind,6]
                                Ls6.set_ydata(datafull[:,6])
                                L6.set_ydata(datafull[:,6])
                            elif Ls7.get_visible():
                                interp = np.linspace(datafull[left_ind,7],datafull[right_ind,7],(right_ind-left_ind))
                                datafull[left_ind:right_ind,7] = datafull[left_ind:right_ind,7] - interp[:] + datafull[left_ind,7]
                                Ls7.set_ydata(datafull[:,7])
                                L7.set_ydata(datafull[:,7])
                            elif Ls8.get_visible():
                                interp = np.linspace(datafull[left_ind,8],datafull[right_ind,8],(right_ind-left_ind))
                                datafull[left_ind:right_ind,8] = datafull[left_ind:right_ind,8] - interp[:] + datafull[left_ind,8]
                                Ls8.set_ydata(datafull[:,8])
                                L8.set_ydata(datafull[:,8])
                            elif Ls9.get_visible():
                                interp = np.linspace(datafull[left_ind,9],datafull[right_ind,9],(right_ind-left_ind))
                                datafull[left_ind:right_ind,9] = datafull[left_ind:right_ind,9] - interp[:] + datafull[left_ind,9]
                                Ls9.set_ydata(datafull[:,9])
                                L9.set_ydata(datafull[:,9])
                            curr_task = 0
    
                        elif curr_task == 1:
                            if x_pos1 > x_pos2:
                                right_ind = int(x_pos1)
                                left_ind = int(x_pos2)
                            else:    
                                right_ind = int(x_pos2)
                                left_ind = int(x_pos1)
                            if Ls1.get_visible():
                                datafull[left_ind:right_ind,1] = datafull[left_ind,1]
                                Ls1.set_ydata(datafull[:,1])
                                L1.set_ydata(datafull[:,1])
                            elif Ls2.get_visible():
                                datafull[left_ind:right_ind,2] = datafull[left_ind,2]
                                Ls2.set_ydata(datafull[:,2])
                                L2.set_ydata(datafull[:,2])
                            elif Ls3.get_visible():
                                datafull[left_ind:right_ind,3] = datafull[left_ind,3]
                                Ls3.set_ydata(datafull[:,3])
                                L3.set_ydata(datafull[:,3])
                            elif Ls4.get_visible():
                                datafull[left_ind:right_ind,4] = datafull[left_ind,4]
                                Ls4.set_ydata(datafull[:,4])
                                L4.set_ydata(datafull[:,4])
                            elif Ls5.get_visible():
                                datafull[left_ind:right_ind,5] = datafull[left_ind,5]
                                Ls5.set_ydata(datafull[:,5])
                                L5.set_ydata(datafull[:,5])
                            elif Ls6.get_visible():
                                datafull[left_ind:right_ind,6] = datafull[left_ind,6]
                                Ls6.set_ydata(datafull[:,6])
                                L6.set_ydata(datafull[:,6])
                            elif Ls7.get_visible():
                                datafull[left_ind:right_ind,7] = datafull[left_ind,7]
                                Ls7.set_ydata(datafull[:,7])
                                L7.set_ydata(datafull[:,7])
                            elif Ls8.get_visible():
                                datafull[left_ind:right_ind,8] = datafull[left_ind,8]
                                Ls8.set_ydata(datafull[:,8])
                                L8.set_ydata(datafull[:,8])
                            elif Ls9.get_visible():
                                datafull[left_ind:right_ind,9] = datafull[left_ind,9]
                                Ls9.set_ydata(datafull[:,9])
                                L9.set_ydata(datafull[:,9])
                            curr_task = 0
                            
                        elif curr_task == 2:
                            if x_pos1 > x_pos2:
                                right_ind = int(x_pos1)
                                left_ind = int(x_pos2)
                            else:    
                                right_ind = int(x_pos2)
                                left_ind = int(x_pos1)
                            if Ls1.get_visible():
                                datafull[left_ind:right_ind,1] = datafull[right_ind,1]
                                Ls1.set_ydata(datafull[:,1])
                                L1.set_ydata(datafull[:,1])
                            elif Ls2.get_visible():
                                datafull[left_ind:right_ind,2] = datafull[right_ind,2]
                                Ls2.set_ydata(datafull[:,2])
                                L2.set_ydata(datafull[:,2])
                            elif Ls3.get_visible():
                                datafull[left_ind:right_ind,3] = datafull[right_ind,3]
                                Ls3.set_ydata(datafull[:,3])
                                L3.set_ydata(datafull[:,3])
                            elif Ls4.get_visible():
                                datafull[left_ind:right_ind,4] = datafull[right_ind,4]
                                Ls4.set_ydata(datafull[:,4])
                                L4.set_ydata(datafull[:,4])
                            elif Ls5.get_visible():
                                datafull[left_ind:right_ind,5] = datafull[right_ind,5]
                                Ls5.set_ydata(datafull[:,5])
                                L5.set_ydata(datafull[:,5])
                            elif Ls6.get_visible():
                                datafull[left_ind:right_ind,6] = datafull[right_ind,6]
                                Ls6.set_ydata(datafull[:,6])
                                L6.set_ydata(datafull[:,6])
                            elif Ls7.get_visible():
                                datafull[left_ind:right_ind,7] = datafull[right_ind,7]
                                Ls7.set_ydata(datafull[:,7])
                                L7.set_ydata(datafull[:,7])
                            elif Ls8.get_visible():
                                datafull[left_ind:right_ind,8] = datafull[right_ind,8]
                                Ls8.set_ydata(datafull[:,8])
                                L8.set_ydata(datafull[:,8])
                            elif Ls9.get_visible():
                                datafull[left_ind:right_ind,9] = datafull[right_ind,9]
                                Ls9.set_ydata(datafull[:,9])
                                L9.set_ydata(datafull[:,9])
                            curr_task = 0
    
                        ax3.figure.canvas.draw() 
    
    
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    cid = fig.canvas.mpl_connect('pick_event', onpick)

def datacorr(filename):

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
    global yvar
    global xvar    
    global txtf1
    global txtf2
    global F_ind
    global is_marked

    data = np.genfromtxt(filename, delimiter=",")
    
    max_slices = 5
    
    cc = np.shape(data)[1]
    cr = np.shape(data)[0]
    
    datafull = np.empty((np.shape(data)[0],10))
    datafull[:,:] = np.NaN
    datafullO = np.empty((np.shape(data)[0],10))
    marked = np.zeros((np.shape(data)[0],1))
    dataOut = np.empty((cr,cc))
    
    fig = plt.figure()
    ax1 = plt.subplot2grid((3,3), (0,0), colspan=1)
    ax2 = plt.subplot2grid((3,3), (0,1), colspan=2)
    ax3 = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan=2)
    
    #ax1 = plt.subplot2grid((2,2), (0,0), colspan=1)
    #ax2 = plt.subplot2grid((2,2), (0,1), colspan=2)
    #ax3 = plt.subplot2grid((2,2), (1,1), colspan=2, rowspan=2)
    
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
    #ax3.axes.set_ylim(-0.1, 1.1)
    
    xvar = 1
    yvar = 2
    
    L1, L2, L3, L4, L5, L6, L7, L8, L9 = plot_lines(ax2, datafull)
    Ls1, = ax3.plot(datafullO[:,xvar-1], datafullO[:,yvar-1], 'ko', markersize = 6)#, picker=5)
    Ls2, = ax3.plot(datafullO[0,xvar-1], datafullO[0,yvar-1], 'bo', markersize = 10)#, picker=5)
    Ls2.set_visible(False)
    #Ls2, = ax3.plot(datafull[:,0], datafullO[:,2], picker=5)
    #Ls3, = ax3.plot(datafull[:,0], datafullO[:,3], picker=5)
    #Ls4, = ax3.plot(datafull[:,0], datafullO[:,4], picker=5)
    #Ls5, = ax3.plot(datafull[:,0], datafullO[:,5], picker=5)
    #Ls6, = ax3.plot(datafull[:,0], datafullO[:,6], picker=5)
    #Ls7, = ax3.plot(datafull[:,0], datafullO[:,7], picker=5)
    #Ls8, = ax3.plot(datafull[:,0], datafullO[:,8], picker=5)
    #Ls9, = ax3.plot(datafull[:,0], datafullO[:,9], picker=5)
    
    Ls0, = ax2.plot([-9999, -9999], [-0.1, 1.5], 'k', linestyle='dashed')
    
    x = [1.8, 1.8, 1.8, 1.8, 1.8, 1.8]
    y = [0.6, 1.6, 2.6, 3.6, 4.6, 5.6]
    
    ax1.plot(x, y, 'bo', picker=5, ms=10)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)
    #ax3.axes.get_xaxis().set_visible(False)
    #ax3.axes.get_yaxis().set_visible(False)
    ax1.axes.set_xlim(0, 2)
    ax1.axes.set_ylim(0, 6.1)
    
    sub_set_ind = 0
    sub_set = 1
    
    ax1.text(0.5, 2.8, 'X Var', fontsize=12)
    ax1.text(0.5, 4.8, 'Y Var', fontsize=12)
    ax1.text(1.1, 1.35, 'Mark', fontsize=12)
    ax1.text(1.1, 0.35, 'Exit', fontsize=12)
    
    filtBx1 = [1.8, 1.8]
    filtBy1 = [3.6, 2.6]
    filtBx2 = [1.8, 1.8]
    filtBy2 = [5.6, 4.6]
    
    txtf1 = ax1.text(1.2, 2.8, xvar, fontsize=12)
    ax1.plot(filtBx1[0], filtBy1[0], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx1[1], filtBy1[1], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx1[0], filtBy1[0], 'b^', ms=8)
    ax1.plot(filtBx1[1], filtBy1[1], 'bv', ms=8)
    
    txtf2 = ax1.text(1.2, 4.8, yvar, fontsize=12)
    ax1.plot(filtBx2[0], filtBy2[0], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx2[1], filtBy2[1], 'wo', markeredgewidth=0, ms=14)
    ax1.plot(filtBx2[0], filtBy2[0], 'b^', ms=8)
    ax1.plot(filtBx2[1], filtBy2[1], 'bv', ms=8)
    
    #ax2.text(np.round(np.shape(data)[0]/2,0), 1.1, 'Full data set', fontsize=12, ha='center')
    ax3.text(np.round(np.shape(data)[0]/2,0), 1.0, 'Correlation', fontsize=12, ha='center')
    
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
    right_pos = int(ind_max)
    first_click = 0
    last_click = 'left'
    zooming = 0
    curr_slice = 1
    is_marked = 0
    
    ax2.axes.set_xlim(0, ind_max)         
    #ax3.margins(x=0.1, y=0.1)
    #ax3.axes.set_xlim(0, ind_max)         
    
       
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
        global yvar
        global xvar    
        global txtf1
        global txtf2
        global F_ind
        global is_marked
        
        mouseevent = event.mouseevent
        axs = mouseevent.inaxes
    #    thisline = event.artist
        ind = event.ind
        if axs == ax1:
    #        print(ind)
            if ind == 0:
    
    #                a = np.array([datafullO[left_pos:right_pos,yvar-1], datafullO[left_pos:right_pos,xvar-1]])
                for i in range(1,cc+1):
                    dataOut[:,i-1] = datafullO[:,i-1]
                newdata = np.concatenate((dataOut, marked), axis=1)
                filenamen = filename[0:filename.find('.')] + '_C' + filename[filename.find('.'):len(filename)]
                np.savetxt(filenamen, newdata, fmt='%1.5f', delimiter=",")
                plt.close()
            elif ind == 1:
                if is_marked == 1:
                    marked[F_ind] = 1
                else:
                    marked[left_pos:right_pos] = 1
            elif ind == 5:
                yvar = yvar + 1
                if yvar == 0:
                    yvar = 1
                elif yvar > cc:
                    yvar = cc
                txtf2.remove()
                txtf2 = ax1.text(1.2, 4.8, yvar, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
                Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                is_marked = 0
                Ls2.set_visible(False)
                xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
            elif ind == 4:
                yvar = yvar - 1
                if yvar == 0:
                    yvar = 1
                elif yvar > cc:
                    yvar = cc
                txtf2.remove()
                txtf2 = ax1.text(1.2, 4.8, yvar, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
                Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                is_marked = 0
                Ls2.set_visible(False)
                xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
            elif ind == 3:
                xvar = xvar + 1
                if xvar == 0:
                    xvar = 1
                elif xvar > cc:
                    xvar = cc
                txtf1.remove()
                txtf1 = ax1.text(1.2, 2.8, xvar, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
                Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                is_marked = 0
                Ls2.set_visible(False)
                xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
            elif ind == 2:
                xvar = xvar - 1
                if xvar == 0:
                    xvar = 1
                elif xvar > cc:
                    xvar = cc
                txtf1.remove()
                txtf1 = ax1.text(1.2, 2.8, xvar, fontsize=12)
                ax1.figure.canvas.draw() # this line is critical to change the linewidth
                Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                is_marked = 0
                Ls2.set_visible(False)
                xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
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
        global yvar
        global xvar    
        global F_ind
        global is_marked
    
        if (fig.canvas.manager.toolbar._active is None) == True:
    
        #    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
        #        event.button, event.x, event.y, event.xdata, event.ydata)
            xpos = event.xdata
            ypos = event.ydata
            axs = event.inaxes
            if axs == ax2:
        #        plt.axvline(x=1) 
        #        xl = [xpos, xpos]
        #        yl = [0, 1000]
        #        ax2.axvline(xpos)
                if first_click == 0:
                    leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                    first_click = 1
                    left_pos = int(xpos)
                    left_exist = 1
                    ax2.figure.canvas.draw() # this line is critical to change the linewidth
                    is_marked = 0
                    Ls2.set_visible(False)
                    Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                    Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                    xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                    ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                    xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                    ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                    ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                    ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                    ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                elif first_click == 1:
                    if xpos < left_pos:
                        xpos = left_pos + 10
                    rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                    first_click = 2
                    right_pos = int(xpos)
                    right_exist = 1
                    Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
                    ax2.figure.canvas.draw() # this line is critical to change the linewidth
                    is_marked = 0
                    Ls2.set_visible(False)
                    Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                    Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                    xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                    ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                    xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                    ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                    ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                    ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                    ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                if first_click == 2:
                    middle = (left_pos + right_pos) / 2
                    if xpos < left_pos:
                         leftgrey.remove()
                         leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                         left_pos = int(xpos)
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
                         ax2.figure.canvas.draw() # this line is critical to change the linewidth
                         is_marked = 0
                         Ls2.set_visible(False)
                         Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                         Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                         xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                         ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                         xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                         ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                         ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                         ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                         ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                    elif xpos > right_pos:
                         rightgrey.remove()
                         rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                         right_pos = int(xpos)
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
                         ax2.figure.canvas.draw() # this line is critical to change the linewidth
                         is_marked = 0
                         Ls2.set_visible(False)
                         Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                         Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                         xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                         ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                         xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                         ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                         ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                         ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                         ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                    elif xpos > middle:
                         rightgrey.remove()
                         rightgrey = ax2.axvspan(xpos, ind_max, facecolor='0.5', alpha=0.5)
                         right_pos = int(xpos)
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
                         ax2.figure.canvas.draw() # this line is critical to change the linewidth
                         is_marked = 0
                         Ls2.set_visible(False)
                         Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                         Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                         xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                         ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                         xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                         ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                         ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                         ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                         ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
                    elif xpos < middle: 
                         leftgrey.remove()
                         leftgrey = ax2.axvspan(1, xpos, facecolor='0.5', alpha=0.5)
                         left_pos = int(xpos)
                         Ls0.set_xdata([np.round((right_pos-left_pos)/2+left_pos,0), np.round((right_pos-left_pos)/2+left_pos,0)])            
                         ax2.figure.canvas.draw() # this line is critical to change the linewidth
                         is_marked = 0
                         Ls2.set_visible(False)
                         Ls1.set_ydata(datafullO[left_pos:right_pos,yvar-1])
                         Ls1.set_xdata(datafullO[left_pos:right_pos,xvar-1])
                         xmax = np.max(datafullO[left_pos:right_pos,xvar-1])
                         ymax = np.max(datafullO[left_pos:right_pos,yvar-1])
                         xmin = np.min(datafullO[left_pos:right_pos,xvar-1])
                         ymin = np.min(datafullO[left_pos:right_pos,yvar-1])
                         ax3.axes.set_xlim(xmin-xmax*0.05, xmax+xmax*0.05)
                         ax3.axes.set_ylim(ymin-ymax*0.05, ymax+ymax*0.05)
                         ax3.figure.canvas.draw() # this line is critical to change the linewidth
        
    #            ax3.axes.set_xlim(left_pos, right_pos)         
        #        ax3.axes.set_ylim(-0.1, 1.1)
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
            elif axs == ax3:
                if is_marked == 0:
                    def find_nearest(a,a0):
    #                    idx = (np.abs(array-value)).argmin()
                        idx = np.sum(np.abs(a-a0),1).argmin()
                        return idx
                    
                    a = np.array([datafullO[left_pos:right_pos,yvar-1], datafullO[left_pos:right_pos,xvar-1]])
                    a = a.T
                    a0 = [ypos, xpos]
                    F_ind = find_nearest(a, a0)                    
                    Ls2.set_ydata(datafullO[F_ind,yvar-1])
                    Ls2.set_xdata(datafullO[F_ind,xvar-1])
                    Ls2.set_visible(True)
                    is_marked = 1
                else:
                    Ls2.set_visible(False)
                    ax3.figure.canvas.draw()
    
                    def find_nearest(a,a0):
    #                    idx = (np.abs(array-value)).argmin()
                        idx = np.sum(np.abs(a-a0),1).argmin()
                        return idx
                    a = np.array([datafullO[left_pos:right_pos,yvar-1], datafullO[left_pos:right_pos,xvar-1]])
                    a = a.T
                    a0 = [ypos, xpos]
                    F_ind = find_nearest(a, a0) + left_pos   
                    Ls2.set_ydata(datafullO[F_ind,yvar-1])
                    Ls2.set_xdata(datafullO[F_ind,xvar-1])
                    Ls2.set_visible(True)
                    is_marked = 1
                ax3.figure.canvas.draw() # this line is critical to change the linewidth
    
            
    cid = fig.canvas.mpl_connect('pick_event', onpick)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    
