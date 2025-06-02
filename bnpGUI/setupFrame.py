"""
Created on Tue Aug  3 11:22:02 2021

@author: graceluo

Construct setup frame for bnp_gui
"""
#!/home/beams/USERBNP/.conda/envs/py36/bin/python

import tkinter as tk
from tkinter import ttk
import os, h5py, time
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import colors, patches
from pvComm import pvComm, pvCommsubclass
from misc import coordinate_transform, checkEntryDigit
from skimage import io
import epics as PV
import cv2
import matplotlib.pyplot as plt
class setupFrame:
    
    def Display2Data(self, Axe, x, y):  #return coordinates of a point (x,y)
        return Axe.transData.inverted().transform(np.array([(x, y)]))[0]

    def choose_folder(self):  #choose user folder
        #if self.pvComm.getDir() is not None:
         #   initialDir = os.path.join(self.pvComm.getDir(), "img.dat")  # go to img.dat folder
          #  initialDir = (
           #     initialDir if os.path.exists(initialDir) else self.pvComm.getDir()  #else: get dir till subdir
            #)
        #else:
         #   initialDir = "/mnt/micdata1/2idd"
        self.h5_folder = tk.filedialog.askdirectory(initialdir=pvCommsubclass().user_di()) #create a dialog box that allows the user to select a directory from their file system.the initial is bnp folder
        #------------------xyl: show cycle/user/folder----------------------------------------
        h5_f_component = self.h5_folder.split('/')
        start_idx = h5_f_component.index('bnp') + 1  #get the index of folder after bnp
        self.toshow_folder = os.path.join(*h5_f_component[start_idx:start_idx+2])  #to show folder path after bnp
        self.selectfolder.set(self.toshow_folder)   #in __init__, selectfolder is string variable, sets to h5_folder now
        ext = (".h5", ".tif", ".tiff")  #xyl: revise to add ptycho images
        if len(self.h5_folder) > 10:
            #self.files_combobox["values"] = [
             #   i for i in os.listdir(self.h5_folder) if i[-len(ext) :] == ext
            #]
            self.file_paths = [f for f in os.listdir(self.h5_folder) if f.endswith(ext)]
            
            # Update the file selection dropdown with the list of files
            self.files_combobox.configure(values=self.file_paths)
            self.file_paths.sort(key=lambda f: os.path.getctime(os.path.join(self.h5_folder, f)),reverse=True)

        else:
            self.h5_folder = tk.filedialog.askdirectory(
                initialdir=pvCommsubclass().user_di()
            )
        self.enter_h5_scan.set(pvCommsubclass().next_scan_num())
        
    def updateFileList(self):
        ext = (".h5", "tif", "tiff") #xyl: update to add ptycho image
        #self.files_combobox["values"] = [
            #i for i in os.listdir(self.h5_folder) if i[-len(ext) :] == ext
        #]
        self.file_paths = [f for f in os.listdir(self.h5_folder) if f.endswith(ext)]
        self.files_combobox.configure(values=self.file_paths)
        self.file_paths.sort(key=lambda f: os.path.getctime(os.path.join(self.h5_folder, f)),reverse=True)
        self.enter_h5_scan.set(pvCommsubclass().next_scan_num())
    '''
    def load_scan(self, *args):
        self.h5_filename = os.path.join(self.h5_folder, self.files_combobox.get())
        try:
            self.h5.close()
        except:
            pass

        try:
            self.ShowROI_Rectangle.remove()
        except (AttributeError, ValueError):
            pass

        t_diff = 10
        fmtime = os.path.getmtime(self.h5_filename)
        ctime = time.time()
        if (ctime - fmtime) >= t_diff:
            try:
                self.h5 = h5py.File(self.h5_filename, "r")
                dets = []
                i_det = max(0, self.detector_combobox.current())
                dets = (
                    self.h5["/MAPS/channel_names"][:].astype(str).tolist()
                    + self.h5["/MAPS/scaler_names"][:].astype(str).tolist()
                )
                elmScalers = np.vstack(
                    (self.h5["/MAPS/XRF_roi"][:], self.h5["/MAPS/scalers"][:])
                )
                self.detector_combobox["values"] = dets
                self.detector_combobox.current(i_det)
                self.x = self.h5["/MAPS/x_axis"][()]
                self.y = self.h5["/MAPS/y_axis"][()]
                pvlist = self.h5["/MAPS/extra_pvs"][0].astype(str).tolist()
                pvval = self.h5["/MAPS/extra_pvs"][1].astype(str).tolist()
                if len(pvlist) < 5:
                    self.file_z = None
                    self.file_theta = None
                
                else:
                    self.file_z = float(
                        pvval[pvlist.index(self.pvComm.pvs["z_value_Act"].pv.pvname)]
                    )
                    self.file_theta = float(
                        pvval[pvlist.index(self.pvComm.pvs["sm_rot_Act"].pv.pvname)]
                    )
                
                self.Image2D = self.Axe2D.imshow(
                    elmScalers[i_det],
                    aspect="equal",
                    interpolation="nearest",
                    cmap="inferno",
                    origin="lower",
                    #                                                 extend = (np.min(self.x),np.max(self.x), np.min(self.y), np.max(self.y)),
                    norm=colors.SymLogNorm(linthresh=0.5)
                    if self.log_button.config("relief")[-1] == "sunken"
                    else colors.Normalize(),
                )
                self.Canvas2D.draw()
                if self.file_z is None:
                    self.openfilemsg.set("Samz PV not found %s" % self.h5_filename)
                    self.open_msg_label.config(fg="red")
                else:
                    self.openfilemsg.set("%s is open" % self.h5_filename)
                    self.open_msg_label.config(fg="green")
            except:
                self.openfilemsg.set("Having trouble opening %s" % self.h5_filename)
                self.open_msg_label.config(fg="red")
        else:
            self.openfilemsg.set("File not ready, getting update from other process")
            self.open_msg_label.config(fg="red")

    def plot_data(self, *args):
        i_det = self.detector_combobox.current()
        elmScalers = np.vstack(
            (self.h5["/MAPS/XRF_roi"][:], self.h5["/MAPS/scalers"][:])
        )
        self.Image2D.set_array(np.array(elmScalers[i_det]))
        if self.log_button.config("relief")[-1] == "sunken":
            # self.Image2D.set_norm(colors.LogNorm())
            self.Image2D.set_norm(colors.SymLogNorm(linthresh=0.5))
        else:
            self.Image2D.set_norm(colors.Normalize())
        self.Canvas2D.draw()
    '''
    def load_scan(self, *args):
        self.h5_filename = os.path.join(self.h5_folder, self.files_combobox.get())   
        nameshow = os.path.basename(self.h5_filename)  #just for message showing file name
        try:
            self.h5.close()
            self.Canvas2D.delete("all")
        except:
            pass
    
        try:
            self.ShowROI_Rectangle.remove()
        except (AttributeError, ValueError):
            pass

        if self.h5_filename.endswith('.h5'):
            t_diff = 10
            fmtime = os.path.getmtime(self.h5_filename)
            ctime = time.time()
            if (ctime - fmtime) >= t_diff:
                try:
                    self.h5 = h5py.File(self.h5_filename, "r")
                    dets = []
                    i_det = max(0, self.detector_combobox.current())  #detector_combobox.current() is to get the index of the selected item
                    dets = (
                        self.h5["/MAPS/channel_names"][:].astype(str).tolist() # channel name list
                        + self.h5["/MAPS/scaler_names"][:].astype(str).tolist() #scaler name list
                    ) 
                    elmScalers = np.vstack(
                        (self.h5["/MAPS/XRF_roi"][:], self.h5["/MAPS/scalers"][:])
                    )
                    self.detector_combobox["values"] = dets
                    self.detector_combobox.current(i_det)
                    self.x = self.h5["/MAPS/x_axis"][()]
                    self.y = self.h5["/MAPS/y_axis"][()]
                    pvlist = self.h5["/MAPS/extra_pvs"][0].astype(str).tolist()
                    pvval = self.h5["/MAPS/extra_pvs"][1].astype(str).tolist()
                    if len(pvlist) < 5:
                        self.file_z = None
                        self.file_theta = None
                    
                    else:
                        self.file_z = float(
                            pvval[pvlist.index(self.pvComm.pvs["z_value_Act"].pv.pvname)]
                        )
                        self.file_theta = float(
                            pvval[pvlist.index(self.pvComm.pvs["sm_rot_Act"].pv.pvname)]
                        )
                
                    self.Image2D = self.Axe2D.imshow(
                        elmScalers[i_det],
                        aspect="equal",
                        interpolation="nearest",
                        cmap="viridis",  #change to another cmap, maybe add cmap dropdow box later
                        origin="lower",
                        #                                                 extend = (np.min(self.x),np.max(self.x), np.min(self.y), np.max(self.y)),
                        norm=colors.SymLogNorm(linthresh=0.5)
                        if self.log_button.config("relief")[-1] == "sunken"
                        else colors.Normalize(),
                    )
                    self.Canvas2D.draw()
                    if self.file_z is None:
                        self.openfilemsg.set("Samz PV not found %s" % nameshow)
                        self.open_msg_label.config(fg="red")
                    else:
                        self.openfilemsg.set("%s is open" % nameshow)
                        self.open_msg_label.config(fg="green")
                except:
                    self.openfilemsg.set("Having trouble opening %s" % nameshow)
                    self.open_msg_label.config(fg="red")
            else:
                self.openfilemsg.set("File not ready, getting update from other process")
                self.open_msg_label.config(fg="red")    
        elif self.h5_filename.endswith('.tif') or self.h5_filename.endswith('.tiff'):
            h5_f = pvCommsubclass().user_h5_folder()
            h5_pty_num = self.enter_h5_scan.get()
            h5_pty_num = str(h5_pty_num).zfill(4)
            h5_pty_sc = f'bnp_fly{h5_pty_num}.mda.h5'
            flag = 1   # only for test, change to 0 for real
            if flag == 0:
                self.simpty_h5 = os.path.join(h5_f, h5_pty_sc)     #get the h5 scan for pty
            else:   #for test
                h5_test = '/mnt/micdata1/bnp/2023-1/Isaure/img.dat'
                self.simpty_h5 = os.path.join(h5_test, h5_pty_sc) 
            with h5py.File(self.simpty_h5, 'r') as dat:    
                self.x = dat["/MAPS/x_axis"][()]
                self.y = dat["/MAPS/y_axis"][()]  #get the stage x and y, but ptycho image shape is different from stage x and y move
                self.file_z = np.float16(np.char.decode(dat['/MAPS/extra_pvs'][:])[1][9])   #z read from h5 file 
                self.file_theta = np.float16(np.char.decode(dat['/MAPS/extra_pvs'][:])[1][8])  #theta from h5
            st_x = self.x[0]
            st_y = self.y[0]
            ed_x = self.x[-1]
            ed_y = self.y[-1]
            self.log_button.config(state="disabled")  #to disable the log check box when ptycho image
            self.ptycho_img = io.imread(self.h5_filename)
            self.ptycho_img_yshape = self.ptycho_img.shape[0]
            self.ptycho_img_xshape = self.ptycho_img.shape[1]
            self.x = np.linspace(st_x, ed_x, self.ptycho_img_xshape)
            self.y = np.linspace(st_y, ed_y, self.ptycho_img_yshape)
             

            
            self.Image2D = self.Axe2D.imshow(self.ptycho_img,aspect="equal",
                interpolation="nearest",
                cmap="viridis",  #change to another cmap, maybe add cmap dropdow box later
                origin="lower",
                #                                                 extend = (np.min(self.x),np.max(self.x), np.min(self.y), np.max(self.y)),
                norm=colors.SymLogNorm(linthresh=0.5)
                if self.log_button.config("relief")[-1] == "sunken"
                else colors.Normalize(),
            )
            self.Canvas2D.draw()
            self.openfilemsg.set("%s is open" % nameshow)
            self.open_msg_label.config(fg="green")
  #-----------------------xyl: for compare pty and xrf images----------------------------------          
    def compare_pty_xrf(self):
        flag = 1   # for test, change to 0 for real
        h5_pty_num = self.enter_h5_scan.get()
        h5_pty_num = str(h5_pty_num).zfill(4)
        xrf_elm = self.elm_pty.get()
        xrf_elm2 = xrf_elm.split(',')
        fig, axes = plt.subplots(2,2)
        fig2,ax11 = plt.subplots(1,1)
        axes[0,0].imshow(self.ptycho_img, cmap='gray',origin='lower')
        axes[0,0].set_title('pty')
        h5_pty_sc = f'bnp_fly{h5_pty_num}.mda.h5'
        if flag == 0:
            h5_f = pvCommsubclass().user_h5_folder()                       
            simpty_h5 = os.path.join(h5_f, h5_pty_sc) 
        else:
            h5_test = '/mnt/micdata1/bnp/2023-1/Isaure/img.dat'   
            simpty_h5 = os.path.join(h5_test, h5_pty_sc) 
        for xe,ax in zip(xrf_elm2,axes.flatten()[1:]):
            with h5py.File(simpty_h5, 'r') as dat:
                xrfdata = dat['/MAPS/XRF_roi'][:]
                channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
                elm_idx = np.where(channel_names==xe)[0][0]
                elmmap = xrfdata[elm_idx,:,:]
                ax.imshow(elmmap,cmap='viridis',origin='lower',vmin=np.min(elmmap),vmax=np.max(elmmap))
                ax.set_title(f'{xe}')
        extent = 0,elmmap.shape[1],0,elmmap.shape[0]
        ax11.imshow(elmmap,cmap='viridis',origin='lower',alpha=0.7,vmin=0,vmax=np.max(elmmap)-1,extent=extent)
        ax11.imshow(self.ptycho_img,cmap='gray',origin='lower',alpha=0.3,vmin=np.min(self.ptycho_img),vmax=np.max(self.ptycho_img),extent=extent)
        ax11.set_title(f'pty+{xrf_elm2[-1]}')
        fig.tight_layout()  
        plt.show()
        
    
    def plot_data(self, *args):
        if self.h5_filename.endswith('.h5'):
            i_det = self.detector_combobox.current()
            elmScalers = np.vstack(
                (self.h5["/MAPS/XRF_roi"][:], self.h5["/MAPS/scalers"][:])
            )
            self.Image2D.set_array(np.array(elmScalers[i_det]))
            if self.log_button.config("relief")[-1] == "sunken":
                # self.Image2D.set_norm(colors.LogNorm())
                self.Image2D.set_norm(colors.SymLogNorm(linthresh=0.5))
            else:
                self.Image2D.set_norm(colors.Normalize())
        elif self.h5_filename.endswith('.tif') or self.h5_filename.endswith('.tiff'):
            self.Image2D.set_array(self.ptycho_img)
        self.Canvas2D.draw()
        
    def logscale_changed(self):
        if self.log_button.config("relief")[-1] == "sunken":
            self.log_button.config(relief="raised")
            self.Image2D.set_norm(colors.Normalize())
        else:
            self.log_button.config(relief="sunken")
            self.Image2D.set_norm(colors.SymLogNorm(linthresh=0.5))
        self.Canvas2D.draw()

    def Draw_Rectangle(self, xmin, xmax, ymin, ymax, color, lw=2, animated=False):
        Rectangle = patches.Rectangle(
            (xmin, ymin),
            width=xmax - xmin,
            height=ymax - ymin,
            alpha=1,
            edgecolor="w",
            fill=False,
            linewidth=lw,
            animated=animated,
        )
        self.Axe2D.add_patch(Rectangle)
        return Rectangle

    '''
    def Canvas2D_Button_Pressed(self, event):
        self.xstart, self.ystart = list(
            map(
                lambda x: int(round(x, 0)),
                (self.Display2Data(self.Axe2D, event.x, event.y)),
            )
        )
        if event.button == 1:
            self.xycorr.set(
                "x, y: (%.2f, %.2f)" % (self.x[self.xstart], self.y[self.ystart])
            )
        if (event.button == 1) & (self.insertType.get() == "XY Center"):
            try:
                self.ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.updateXYcenter()
        elif (event.button == 3) & (self.insertType.get() == "ScanBox"):
            self.Rectangle_Drawing = True
            try:
                self.ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.Canvas2D.mpl_disconnect(self.Canvas2D_Button_Press_Event)
            self.Canvas2D_Button_Release_Event = self.Canvas2D.mpl_connect(
                "button_release_event", self.Canvas2D_Button_Released
            )
     '''  
    #---------------------------xyl: add event for ptycho image--------------------------------        
    def Canvas2D_Button_Pressed(self, event):        
        self.xstart, self.ystart = list(
            map(
                lambda x: int(round(x, 0)),
                (self.Display2Data(self.Axe2D, event.x, event.y)),
            ))  #to get the pixel coordinates of the start point
        
        if event.button == 1:
            self.xycorr.set(
                "x, y: (%.2f, %.2f)" % (self.x[self.xstart], self.y[self.ystart])
            )
        if (event.button == 1) & (self.insertType.get() == "XY Center"):
            try:
                self.ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.updateXYcenter()
        
        elif (event.button == 3) & (self.insertType.get() == "ScanBox"):
            self.Rectangle_Drawing = True
            try:
                self.ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            self.Canvas2D.mpl_disconnect(self.Canvas2D_Button_Press_Event)
            self.Canvas2D_Button_Release_Event = self.Canvas2D.mpl_connect(
                "button_release_event", self.Canvas2D_Button_Released
            ) 
        
          
    def Canvas2D_Button_Released(self, event):
        self.xend, self.yend = list(
            map(
                lambda x: int(round(x, 0)),
                (self.Display2Data(self.Axe2D, event.x, event.y)),
            )
        )

        self.Canvas2D.mpl_disconnect(self.Canvas2D_Button_Release_Event)
        self.Canvas2D_Button_Press_Event = self.Canvas2D.mpl_connect(
            "button_press_event", self.Canvas2D_Button_Pressed
        )

        if self.Rectangle_Drawing:
            self.Rectangle_Drawing = False
            self.non_animated_background = None
            try:
                self.ShowROI_Rectangle.remove()
            except (AttributeError, ValueError):
                pass
            xmin = int(min(self.xend, self.xstart))
            xmax = int(max(self.xend, self.xstart))
            ymin = int(min(self.yend, self.ystart))
            ymax = int(max(self.yend, self.ystart))
            self.ShowROI_Rectangle = self.Draw_Rectangle(xmin, xmax, ymin, ymax, "w")
            self.updateScanBoxParms()
    '''
    def Canvas2D_Mouser_Hover(self, event):
        self.xend, self.yend = list(
            map(
                lambda x: int(round(x, 0)),
                (self.Display2Data(self.Axe2D, event.x, event.y)),
            )
        )
        try:
            dim_y, dim_x = self.Image2D.get_array().shape
        except:
            return
        if (
            (self.yend < dim_y - 0.5)
            and (self.yend > -0.5)
            and (self.xend > -0.5)
            and (self.xend < dim_x - 0.5)
        ):
            if self.Rectangle_Drawing:
                if self.non_animated_background != None:
                    # restore the clean slate background
                    self.Canvas2D.restore_region(self.non_animated_background)
                    if self.xstart > self.xend:  # modify the starting point
                        self.Rectangle.set_x(self.xend)
                    self.Rectangle.set_width(abs(self.xend - self.xstart))
                    if self.ystart > self.yend:  # modify the starting point
                        self.Rectangle.set_y(self.yend)
                    self.Rectangle.set_height(abs(self.yend - self.ystart))
                    self.Axe2D.draw_artist(self.Rectangle)
                    self.Canvas2D.blit(self.Axe2D.bbox)
                else:
                    xmax = int(max(self.xend, self.xstart))
                    xmin = int(min(self.xend, self.xstart))
                    ymax = int(max(self.yend, self.ystart))
                    ymin = int(min(self.yend, self.ystart))
                    self.Rectangle = self.Draw_Rectangle(
                        xmin, xmax, ymin, ymax, "w", animated=True
                    )
                    self.Canvas2D.draw()
                    self.non_animated_background = self.Canvas2D.copy_from_bbox(
                        self.Axe2D.bbox
                    )
    '''
    def Canvas2D_Mouser_Hover(self, event):
        self.message = ""
        self.xend, self.yend = list(
            map(
                lambda x: int(round(x, 0)),
                (self.Display2Data(self.Axe2D, event.x, event.y)),
            )
        )
        try:
            dim_y, dim_x = self.Image2D.get_array().shape
        except:
            return
        if (
            (self.yend < dim_y - 0.5)
            and (self.yend > -0.5)
            and (self.xend > -0.5)
            and (self.xend < dim_x - 0.5)
        ):
            if self.Rectangle_Drawing:
                if self.non_animated_background != None:
                    # restore the clean slate background
                    self.Canvas2D.restore_region(self.non_animated_background)
                    if self.xstart > self.xend:  # modify the starting point
                        self.Rectangle.set_x(self.xend)
                    self.Rectangle.set_width(abs(self.xend - self.xstart))
                    if self.ystart > self.yend:  # modify the starting point
                        self.Rectangle.set_y(self.yend)
                    self.Rectangle.set_height(abs(self.yend - self.ystart))
                    self.Axe2D.draw_artist(self.Rectangle)
                    self.Canvas2D.blit(self.Axe2D.bbox)
                else:
                    xmax = int(max(self.xend, self.xstart))
                    xmin = int(min(self.xend, self.xstart))
                    ymax = int(max(self.yend, self.ystart))
                    ymin = int(min(self.yend, self.ystart))
                    self.Rectangle = self.Draw_Rectangle(
                        xmin, xmax, ymin, ymax, "w", animated=True
                    )
                    self.Canvas2D.draw()
                    self.non_animated_background = self.Canvas2D.copy_from_bbox(
                        self.Axe2D.bbox
                    )
               
    def updateXYcenter(self):
        x_scan = np.round(self.x[self.xstart], 2)
        y_scan = np.round(self.y[self.ystart], 2)
        z_scan = np.round(self.file_z, 2)
        target_theta = self.file_theta

        slabel = ["x_scan", "y_scan", "z_scan", "target_theta"]
        for s in slabel:
            self.scanParms[s].delete(0, tk.END)
            self.scanParms[s].insert(0, "%.2f" % (eval(s)))

        slabel = ["x_theta0", "y_theta0", "z_theta0"]
        if (self.file_theta ** 2) < 1e-3:
            svlabel = ["x_scan", "y_scan", "z_scan"]
        else:
            e = ""
            svlabel = ["e"] * 3

        for s_, sv_ in zip(slabel, svlabel):
            self.scanParms[s_].delete(0, tk.END)
            self.scanParms[s_].insert(0, "%s" % (str(eval(sv_))))

    def updateScanBoxParms(self):
        x_scan = np.round((self.x[self.xstart] + self.x[self.xend]) / 2, 2)
        y_scan = np.round((self.y[self.ystart] + self.y[self.yend]) / 2, 2)
        width = abs(self.x[self.xstart] - self.x[self.xend])
        height = abs(self.y[self.ystart] - self.y[self.yend])
        if self.file_z is not None:    
            z_scan = np.round(self.file_z, 2)
            target_theta = self.file_theta
            slabel = ["x_scan", "y_scan", "width", "height", "z_scan", "target_theta"]
            for s in slabel:
                self.scanParms[s].delete(0, tk.END)
                self.scanParms[s].insert(0, "%.2f" % (eval(s)))
    
            slabel = ["x_theta0", "y_theta0", "z_theta0"]
            if (self.file_theta ** 2) < 1e-3:
                svlabel = ["x_scan", "y_scan", "z_scan"]
            else:
                e = ""
                svlabel = ["e"] * 3
    
            for s_, sv_ in zip(slabel, svlabel):
                self.scanParms[s_].delete(0, tk.END)
                self.scanParms[s_].insert(0, "%s" % (str(eval(sv_))))
        else:
            slabel = ["x_scan", "y_scan", "width", "height"]
            for s in slabel:
                self.scanParms[s].delete(0, tk.END)
                self.scanParms[s].insert(0, "%.2f" % (eval(s)))
            flabel = ['z_scan', 'target_theta', 'x_theta0', 'y_theta0', 'z_theta0']
            for f in flabel: 
                self.scanParms[f].delete(0, tk.END)

    # def checkEntryDigit(self, P):
    #     if (P == "") | (P == "-"):
    #         return True
    #     try:
    #         float(P)
    #         return True
    #     except:
    #         return False
    
    def calcTime(self):
        scant = self.scanType.get()
        if scant == "XANES (fixed region)":
            strlist = self.inputs_labels[2]
            try:
                value = [float(self.scanParms[s].get()) for s in strlist]            
            except:
                value = [0] * len(strlist)
            if 0 in value:
                self.calctime.set("wrong parameters")
                self.calctime_out.config(fg="red")
            else:
                eta_ms = value[3] * value[2] / value[1]
                eta_min = eta_ms / 60 / 0.8
                self.calctime.set("%.3f" % (eta_min))
                self.calctime_out.config(fg="green")
        else:    
            strlist = ["width", "height", "w_step", "h_step", "dwell"]
            try:
                value = [float(self.scanParms[s].get()) for s in strlist]
            except:
                value = [0] * len(strlist)
        
            if value[0] > 80:
                self.calctime.set("Width can not be bigger than 80")
                self.calctime_out.config(fg="red")
            elif 0 not in value:
                eta_ms = value[-1] * value[0] * value[1] / value[2] / value[3]
                eta_min = eta_ms / 1e3 / 60 / 0.8
                self.calctime.set("%.3f" % (eta_min))
                self.calctime_out.config(fg="green")
    #-------------------cal tomo time---------------------------
    def cal_tomo_time(self):
        strlist_tomo = ["width", "height", "w_step", "h_step", "dwell", #0,1,2,3,4
                      "theta_min", "theta_max","theta_inc", #5,6,7
                      "width_fine", "w_step_fine", "height_fine","h_step_fine","dwell_fine"]  #8,9,10,11,12
        try:
            value = [float(self.scanParms[s].get()) for s in strlist_tomo]
        except:
            value = [0] * len(strlist_tomo)
        
        if value[0] > 80 or value[8] > 80 or value[5] <-90 or value[6] > 90:
            self.tomo_est_time.config(text="Width <80, theta=[-90,90]", fg="red")
        elif value[2] == 0 or value[3] == 0 or value[9] == 0 or value[11] == 0 or value[7] == 0:
            self.tomo_est_time.config(text="step size or angle step is 0", fg="red")
        else:
            start = round(value[5])
            end = round(value[6])
            inc_ang = round(value[7])
            tot_ang = round((end - start) / inc_ang) + 1
            one_cs_tm = value[4] * value[0] * value[1] / value[2] / value[3] / 0.8 # dwell * width * height / w_step / h_step in ms unit for one coarse scan overhead 20%
            one_fn_tm = value[12] * value[8] * value[10] / value[9] / value[11] / 0.8 #dwell_fine * width_fine * height_fine / w_step_fine / h_step_fine in ms unit for one coarse scan overhead 20%
            tomo_tm = (one_cs_tm + one_fn_tm) * tot_ang / 1e3 / 60 / 60  # in hrs
            tomo_tm = np.round(tomo_tm,5)
            self.tomo_est_time.config(text=f'Total tomo time (hrs): {tomo_tm}',fg="green")  
   #----------------------------xyl: modify to add time estimation on tomo----------------------------             
    ''' 
   def calcTime(self):
        
        strlist = ["width", "height", "w_step", "h_step", "dwell",
                      "theta_min", "theta_max","theta_inc",
                      "width_fine", "w_step_fine", "height_fine","h_step_fine","dwell_fine"]
        try:
            value = [float(self.scanParms[s].get()) for s in strlist]
        except:
            value = [0] * len(strlist)
    
        if value[0] > 80 or value[8] > 80 or value[5] <-90 or value[6] > 90:
            self.calctime.set("Width <80, theta=[-90,90]")
            self.calctime_out.config(fg="red")
        elif 0 not in value:
            if self.scanType != 'Coarse-Fine':
                eta_ms = value[-1] * value[0] * value[1] / value[2] / value[3]
                eta_min = eta_ms / 1e3 / 60 / 0.8
                self.calctime.set("%.3f" % (eta_min))
                self.calctime_out.config(fg="green")
            else:
                eta_ms = value[-1] * value[0] * value[1] / value[2] / value[3]  #one coarse scan
                eta_min = eta_ms / 1e3 / 60 / 0.8 #20% overhead, for one coarse scan   
                
        else:
            self.calctime.set("setup wrong")
            self.calctime_out.config(fg="red")
    '''

    def getUserDir(self):
        self.pardir.set(self.pvComm.getDir())

    def addToScanBtn(self, func):
        add_scan_btn = tk.Button(
            self.setupfrm, text="Add to scan list", command=func, width=50
        )
        add_scan_btn.grid(row=self.row + 5, column=self.col, columnspan=4)

    def updatexyz0(self):
        if (self.pvComm.getSMAngle() ** 2) < 1e-3:
            self.fillxyz0(empty=False)
            self.coordtform()
            self.updatexyz_msg.set(
                "x-, y-, z- theta0 updated"
                + "\n"
                + "Coordinate transformed. x-, y-, z- scan values updated"
            )
            self.updatexyz_label.config(fg="green")
        else:
            self.updatexyz_msg.set("Theta (Rotation) motor is not at 0")
            self.updatexyz_label.config(fg="red")
            
    def fillxyz0(self, empty = True):
        if empty:
            x,y,z = ['','','']
        else:
            x, y, z = self.pvComm.getXYZcenter()
        for a_, v_ in zip(["x_theta0", "y_theta0", "z_theta0"], [x, y, z]):
            self.scanParms[a_].delete(0, tk.END)
            self.scanParms[a_].insert(0, str(v_)) 
            
    def fillxyzScan(self, empty = True):
        if empty:
            x,y,z,theta = ['','','', '']
        else:
            x, y, z = self.pvComm.getXYZcenter()
            theta = self.pvComm.getSMAngle()
        for a_, v_ in zip(['x_scan', 'y_scan', 'z_scan', 'target_theta'], [x, y, z, theta]):
            self.scanParms[a_].delete(0, tk.END)
            self.scanParms[a_].insert(0, str(v_)) 
    
    def updatexyz(self):
        if (self.pvComm.getSMAngle() ** 2) < 1e-3:
            self.fillxyz0(empty=False)
            self.fillxyzScan()
        else:
            self.fillxyz0()
            self.fillxyzScan(empty=False)


    def coordtform(self):
        fields = ["target_theta", "x_theta0", "y_theta0", "z_theta0"]
        fvals = [float(self.scanParms[a_].get()) for a_ in fields]
        ctform = coordinate_transform(*fvals)
        writeFields = ["x_scan", "y_scan", "z_scan"]
        wfv = [round(ctform[a], 2) for a in ["x", "y", "z"]]
        for a_, v_ in zip(writeFields, wfv):
            self.scanParms[a_].delete(0, tk.END)
            self.scanParms[a_].insert(0, str(v_))
        self.updatexyz_msg.set("Coordinate transformed. x-, y-, z- scan values updated")
        self.updatexyz_label.config(fg="green")
    #------------------------------for export scan images------------------------------------
    def export_scan_images(self):
        flag = 1  #just for test, will remove for real experiment
        if flag ==0:
            value = self.enter_value.get()
            v = value.split(',')
            scan_list = []
            elm = self.enter_value_elm.get()
            savefol = pvCommsubclass().user__fine_folder()
            error_list = []
            if len(v) == 2:
                st = int(v[0])
                ed = int(v[1])
                scan_list = np.arange(st,ed+1,1)
            elif len(v)==3 and v[2] == '2':
                st = int(v[0])
                ed = int(v[1])
                scan_list = np.arange(st,ed+1,2)            
            else:
                scan_list = [int(vv) for vv in v]
            for sl in scan_list:
                sl1 = str(sl).zfill(4)
                file_fine = f'{pvCommsubclass().user_di()}/img.dat/bnp_fly{sl1}.mda.h5'
                try:
                    with h5py.File(file_fine, 'r') as dat:
                        xrfdata = dat['/MAPS/XRF_roi'][:]
                        channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
                        elm_idx = np.where(channel_names==elm)[0][0]   
                        elmmap = xrfdata[elm_idx,:,:]
                        theta = round(float(dat['MAPS/extra_pvs'][1,8].decode('utf-8'))) #this need to be changed based on cycle data
                        theta_fake = theta + 90
                        io.imsave(f'{savefol}/icm{theta_fake}_bnp_fly{sl1}_Ang{theta}_{elm}.tif',elmmap)
                except Exception:
                    error_list.append(sl)
                time.sleep(0.1)
            lenl = len(scan_list)
            if len(error_list) == 0:
                self.export_prog.config(text=f'Total {lenl} scans are exported',fg='green')
            else:
                self.export_prog.config(text=f'Scans {error_list} have issues',fg='red')
        else:
            folder = '/mnt/micdata1/bnp/2022-3/Merk/img.dat.5det.ele'
            value = self.enter_value.get()
            v = value.split(',')
            elm = self.enter_value_elm.get()
            savefol = pvCommsubclass().user__fine_folder()
            error_list = []
            if len(v) == 2:
                st = int(v[0])
                ed = int(v[1])
                scan_list = np.arange(st,ed+1,1)
            elif len(v)==3 and v[2] == '2':
                st = int(v[0])
                ed = int(v[1])
                scan_list = np.arange(st,ed+1,2)
            else:
                scan_list = [int(vv) for vv in v]
            for sl in scan_list:
                sl1 = str(sl).zfill(4)
                file_fine = f'{folder}/bnp_fly{sl1}.mda.h5'                
                try:
                    dat = h5py.File(file_fine,'r')
                    xrfdata = dat['/MAPS/XRF_roi'][:]
                    channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
                    elm_idx = np.where(channel_names==elm)[0][0]   
                    elmmap = xrfdata[elm_idx,:,:]
                    theta = round(float(dat['MAPS/extra_pvs'][1,3].decode('utf-8'))) #this need to be changed based on cycle data
                    theta_fake = theta + 90
                    #io.imshow(elmmap)
                    cv2.imwrite(f'{savefol}/icm{theta_fake}_bnp_fly{sl1}_Ang{theta}_{elm}.tif',elmmap)
                    dat.close()
                    '''
                        with h5py.File(file_fine, 'r') as dat:
                            xrfdata = dat['/MAPS/XRF_roi'][:]
                            channel_names = np.char.decode(dat['MAPS']['channel_names'][:])
                            elm_idx = np.where(channel_names==elm)[0][0]   
                            elmmap = xrfdata[elm_idx,:,:]
                            theta = round(float(dat['MAPS/extra_pvs'][1,3].decode('utf-8'))) #this need to be changed based on cycle data
                            theta_fake = theta + 90
                            io.imsave(f'{savefol}/icm{theta_fake}_bnp_fly{sl1}_Ang{theta}_{elm}.tif',elmmap)
                    '''
                except FileNotFoundError or IndexError:
                    error_list.append(sl)
                    dat.close()
                time.sleep(0.1)
            lenl = len(scan_list)
            if len(error_list) == 0:
                self.export_prog.config(text=f'Total {lenl} scans are exported',fg='green')
            else:
                self.export_prog.config(text=f'Scans {error_list} have issues',fg='red')
                
    #------------------------------create folder------------------------
    def CF_folder(self):
        self.rootfolder = PV.caget('9idbBNP:saveData_fileSystem')  #root till cycle
        self.rootfolder = self.rootfolder.replace('//micdata/data1','/mnt/micdata1')
        self.user = PV.caget('9idbBNP:saveData_subDir').split('/')[0]  #user and mda
        self.user_folder = os.path.join(self.rootfolder,self.user)
        F_dir = 'Fine_images'
        C_dir = 'Coarse_images'
        self.F_folder = os.path.join(self.user_folder,F_dir)
        self.C_folder = os.path.join(self.user_folder,C_dir)        
        if os.path.exists(self.F_folder):
            print('folder to save fine images exists')
        else:
            os.makedirs(self.F_folder)
        if os.path.exists(self.C_folder):
            print('folder to save coarse images exists')
        else:
            os.makedirs(self.C_folder)
        
 

          
    def __init__(self, tabControl):
        self.setupfrm = ttk.Frame(tabControl)  #create a frame inside tabControl

        self.selectfolder = tk.StringVar()  #a mutable string variable
        self.selectfolder.set(" ")  #give the string value to above
        mdafolder_button = tk.Button(
            self.setupfrm, text="Open Folder", command=self.choose_folder
        )  #create the button 'Open Folder' in 'setupfrm' Frame, the function is choose_folder
        mdafolder_button.grid(column=0, row=0, padx=(5, 5), pady=(5, 5))  #setupfrm is in grid layout, mdafolder_button is in 1 col and 1 row
        folder_txt = tk.Label(self.setupfrm, textvariable=self.selectfolder) # create folder_txt label, it display the value of self.selectfolder string
        folder_txt.grid(row=1, column=0,sticky='ew') #size and pos of folder_txt

        #scannum_txt = tk.Label(self.setupfrm, text="Files:")   #create a scannum_txt column, name:'Files:'
        #scannum_txt.grid(column=0, row=53, pady=(5, 5), padx=(5, 5))

        self.scfile_sv = tk.StringVar()
        self.files_combobox = ttk.Combobox(
            self.setupfrm, textvariable=self.scfile_sv, state="readonly", width=20
        )  #a drupdown box, value of scfile_sv is the selection, restrict to the options, cannot input
        self.files_combobox.grid(row=1,column=1, pady=(5, 5))
        self.files_combobox.bind("<<ComboboxSelected>>", self.load_scan)  #<<ComboboxSelected>>" event is generated when the user selects an item from the list
        self.pvComm = pvComm()   #all the pvs
        #        self.scannum_entry = tk.Entry(self.setupfrm, width=10)
        #        self.scannum_entry.grid(column=1, row = 1, padx=(5,5), pady=(5,5))

        loadscan_button = tk.Button(
            self.setupfrm, text="Update", command=self.updateFileList   #update files in folder
        )
        loadscan_button.grid(row=1,column=2, pady=(5, 5),sticky='ew')
        detector_sv = tk.StringVar()
        self.detector_combobox = ttk.Combobox(
            self.setupfrm, textvariable=detector_sv, state="readonly", width=20
        )
        self.detector_combobox.grid(row=1,column=3, pady=(5, 5))
        self.detector_combobox.bind("<<ComboboxSelected>>", self.plot_data)

        self.log_button = tk.Button(
            self.setupfrm, text="Log", command=self.logscale_changed
        )
        self.log_button.grid(column=0, row=52, padx=(4, 5), pady=(0, 0))
        self.xycorr = tk.StringVar()
        self.xycorr.set("x, y: (0.00, 0.00)")
        xycorr_label = tk.Label(self.setupfrm, textvariable=self.xycorr)
        xycorr_label.grid(row=52, column=1)

        self.openfilemsg = tk.StringVar()
        self.openfilemsg.set("")
        self.open_msg_label = tk.Label(self.setupfrm, textvariable=self.openfilemsg)
        self.open_msg_label.grid(row=53, column=0, columnspan=12)
        self.file_theta = 0
        self.file_z = 0

        self.Figure2D = Figure()
        self.Axe2D = self.Figure2D.add_axes([0, 0, 1, 1]) #[0, 0, 1, 1] sets the left and bottom edges of the axes to be at 0% of the figure width and height, respectively, and the width and height of the axes to be 100% of the figure width and height, 
        self.Axe2D.set_axis_off()
        #self.Axe2D.set_axis_on()

        self.Canvas2D = FigureCanvasTkAgg(self.Figure2D, master=self.setupfrm)
        self.Canvas2D.draw()
        self.Canvas2D.get_tk_widget().config(width=600, height=600, cursor="cross")
        self.Canvas2D.get_tk_widget().grid(
            column=0,
            columnspan=10,
            row=2,
            rowspan=50,
            padx=(10, 5),
            pady=(5, 5),
            sticky="W",
        )  #display 2D plot
        #------------------------------------add toolbar-------------------------------------
        toolbar_frame = tk.Frame(master=self.setupfrm)
        toolbar_frame.grid(column=1, row=0, sticky="ew")

        toolbar = NavigationToolbar2Tk(canvas=self.Canvas2D, window=toolbar_frame)
        toolbar.update()
        toolbar.mode = 'pan'
        toolbar.message = ''
        #-----------------------------------add a box to select h5 scan for ptycho image----------------
        self.enter_h5_scan = tk.StringVar()
        h5_scan_for_pty = tk.Label(self.setupfrm, text='h5 scan for pty:',width=13)
        h5_scan_for_pty.grid(row=0,column=2, padx=1,sticky='w')
        h5_scan_for_pty_num = tk.Entry(
                    self.setupfrm,
                    width=5,
                    textvariable=self.enter_h5_scan)
        h5_scan_for_pty_num.grid(row=0,column=3,sticky='w')
        compare_pty_xrf = tk.Button(self.setupfrm, text='Compare pty & xrf',width=13,command=self.compare_pty_xrf)
        compare_pty_xrf.grid(row=52,column=2, sticky='w')
        self.elm_pty = tk.StringVar()
        elm_pty = tk.Entry(
                    self.setupfrm,
                    width=6,
                    textvariable=self.elm_pty)
        elm_pty.grid(row=52,column=3,sticky='w')
        #---------------------------------------------------------------------------------------
        self.Canvas2D_Mouse_Hover_Event = self.Canvas2D.mpl_connect(
            "motion_notify_event", self.Canvas2D_Mouser_Hover
        )
        self.Canvas2D_Button_Press_Event = self.Canvas2D.mpl_connect(
            "button_press_event", self.Canvas2D_Button_Pressed
        )
        self.Rectangle_Drawing = False
        self.Dot_Drawing = False
        self.non_animated_background = None

        self.row = 2
        self.col = 10
        self.scantype_txt = tk.Label(
            self.setupfrm, text="1. Choose a type of scan below:"
        )
        self.scantype_txt.grid(
            row=self.row, column=self.col, padx=(0, 400), columnspan=10
        )

        self.scanType = tk.StringVar(self.setupfrm)
        self.scanType.set("XRF")
        #scantype = ["XRF",'Coarse-Fine (Fixed Angle)', "Angle Sweep", "Coarse-Fine"]
        #--------------------------add XANES-----------------------------------------
        scantype = ["XRF",'Coarse-Fine (Fixed Angle)', "Angle Sweep", "Coarse-Fine","XANES (fixed region)"]
        padx = [(18, 10), (0, 10), (0,10), (18, 10), (0, 10)]
        self.row += 1
        for i, s in enumerate(scantype):
            #----------------------------when click Coarse-Fine, to create a folder to save fine scan image----------
            if s != 'Coarse-Fine':
                a = ttk.Radiobutton(
                    master=self.setupfrm, text=s, 
                    variable=self.scanType, value=s)
            else: 
                a = ttk.Radiobutton(
                    master=self.setupfrm, text=s, 
                    variable=self.scanType, value=s,command=self.CF_folder)                
            #------------------------------------------------------------------------------------------------------------          
            a.grid(row=self.row if i < 3 else (self.row+1), 
                   column=self.col + i%3, padx=padx[i], 
                   stick='w')
        
        # add ptycho checkboc
        self.ptychoVal = tk.IntVar() #interger value
        self.ptychoVal.set(0)
        ptycho_btn = ttk.Checkbutton(
            master=self.setupfrm, text='Ptycho Enabled', variable=self.ptychoVal)
        ptycho_btn.grid(row=self.row+1, column=self.col+2)
        
        self.row += 2
        self.inserttype_txt = tk.Label(
            self.setupfrm, text="2. Choose a method to insert scan parameters:"
        )
        self.inserttype_txt.grid(
            row=self.row, column=self.col, padx=(0, 320), columnspan=8
        )
        self.insertType = tk.StringVar(self.setupfrm)
        self.insertType.set("Manual")
        inserttype = ["Manual", "XY Center", "ScanBox"]
        padx = [(0, 84), (0, 120), (0, 120)]
        self.row += 1
        for i, s in enumerate(inserttype):
            a = ttk.Radiobutton(
                master=self.setupfrm, text=s, variable=self.insertType, value=s
            )
            a.grid(row=self.row, column=self.col + i, padx=padx[i])

        self.scanParms = {}
        self.row += 2
        self.col = 10
        ncol = 3
        '''
        self.inputs_labels = [
            [
                "x_theta0",
                "y_theta0",
                "z_theta0",
                "x_scan",
                "y_scan",
                "z_scan",
                "dwell",
                "width",
                "w_step",
                "target_theta",
                "height",
                "h_step",
            ],
            [
                "theta_min",
                "theta_max",
                "theta_inc",
                "elm",
                "width_fine",
                "w_step_fine",
                "n_clusters",
                "height_fine",
                "h_step_fine",
                "sel_cluster",
                "dwell_fine",
                "use_mask",
                "mask_elm",
            ],
        ]
        '''
        self.inputs_labels = [
            [
                "x_theta0",
                "y_theta0",
                "z_theta0",
                "x_scan",
                "y_scan",
                "z_scan",
                "dwell",
                "width",
                "w_step",
                "target_theta",
                "height",
                "h_step",
            ],
            [
                "theta_min",
                "theta_max",
                "theta_inc",
                "elm",
                "width_fine",
                "w_step_fine",
                "n_clusters",
                "height_fine",
                "h_step_fine",
                "sel_cluster",
                "dwell_fine",
                "use_mask",
                "mask_elm",
                "Gaussian_blur",
                "log_his",
                "flag",   #flag is just for test/real scan change: 0 is real, 1 is test
            ],
            [
                "Energy (keV)",
                "Energy step (keV)",
                "Energy width (keV)",
                "Dwell (s)",
            ]
        ]  #inputs_labels[0]: scan settings on top, inputs_labels[1]: scan settings for coarse-fine, inputs_labels[2]: scan settings for XANES
        #        corr = self.pvComm.getXYZcenter()
        temp_input = [0] * 6 + [20, 20, 1, 15, 10, 1]
        vcmd = self.setupfrm.register(checkEntryDigit)
        for i, s in enumerate(self.inputs_labels[0]):
            if i % ncol == 0:
                self.row += 1
            sv = tk.StringVar()
            sv.trace("w", lambda name, index, mode, sv=sv: self.calcTime())
            tk.Label(self.setupfrm, text=s + ":").grid(
                row=self.row, column=self.col + i % ncol, sticky="W", columnspan=2
            )
            a = tk.Entry(
                self.setupfrm,
                width=10,
                textvariable=sv,
                validate="all",
                validatecommand=(vcmd, "%P"),
            )
            a.grid(row=self.row, column=self.col + i % ncol, padx=(70, 0))
            a.insert(0, str(temp_input[i]))
            #            a.insert(0,'0')
            self.scanParms.update({s: a})

        self.row += 1
        smp_txt = tk.Label(self.setupfrm, text="Sample name:")
        smp_txt.grid(row=self.row, column=self.col, sticky="w")
        self.smp_name = tk.StringVar()
        smp_entry = tk.Entry(self.setupfrm, textvariable=self.smp_name)
        smp_entry.grid(
            row=self.row, column=self.col, columnspan=2, sticky="w", padx=(88, 0)
        )

        bda_txt = tk.Label(self.setupfrm, text="BDA Position:")
        bda_txt.grid(row=self.row, column=self.col + 1, sticky="w", padx=(100, 0))
        self.bda = tk.StringVar()
        bda_entry = tk.Entry(
            self.setupfrm,
            textvariable=self.bda,
            validate="all",
            validatecommand=(vcmd, "%P"),
        )

        # TO DO: get bda position from save state instead from motor
        self.bda.set("%.2f" % (self.pvComm.getBDAx()))
        bda_entry.grid(
            row=self.row, column=self.col + 2, sticky="w", padx=(25, 0), columnspan=2
        )

        self.row += 1
        updatexyz0_btn = tk.Button(
            self.setupfrm, text="XYZ_0 + transform", command=self.updatexyz0
        )
        updatexyz0_btn.grid(row=self.row, column=self.col, columnspan=2)
        updatexyz_btn = tk.Button(
            self.setupfrm, text="XYZ from PV", command=self.updatexyz
        )
        updatexyz_btn.grid(row=self.row, column=self.col + 1, columnspan=2)

        self.row += 1
        calctime_txt = tk.Label(self.setupfrm, text="Estimated scan time (min):")
        calctime_txt.grid(row=self.row, column=self.col, sticky="w")
        self.calctime = tk.StringVar()
        self.calctime.set("0")
        self.calctime_out = tk.Label(
            self.setupfrm, width=18, textvariable=self.calctime
        )
        self.calctime_out.grid(row=self.row, column=self.col + 1, sticky="w")

        self.row += 1
        self.updatexyz_msg = tk.StringVar()
        self.updatexyz_msg.set("")
        self.updatexyz_label = tk.Label(self.setupfrm, textvariable=self.updatexyz_msg)
        self.updatexyz_label.grid(
            row=self.row, column=self.col, sticky="w", columnspan=3
        )

        self.row += 4
        self.inserttype_txt = tk.Label(
            self.setupfrm, text="Additional scan parameters for tomography"
        )
        self.inserttype_txt.grid(
            row=self.row, column=self.col, padx=(0, 0), columnspan=8
        )
        for i, s in enumerate(self.inputs_labels[1]):
            if i % ncol == 0:
                self.row += 1
                    
            tk.Label(self.setupfrm, text=s + ":").grid(
                row=self.row, column=self.col + i % ncol, sticky="W", columnspan=2
            )
            if (s == "elm") | (s == "mask_elm"):
                a = tk.Entry(self.setupfrm, width=10)
                a.insert(0, "P")
            else:
                a = tk.Entry(
                    self.setupfrm,
                    width=10,
                    validate="all",
                    validatecommand=(vcmd, "%P"),
                )
                a.insert(0, "0")
            a.grid(row=self.row, column=self.col + i % ncol, padx=(70, 0))
            self.scanParms.update({s: a})
        # self.scanParms.update({'ptycho':self.ptychoVal})
        

        row = self.row + 7
        #-------------------------add a label/button for tomo time estimation------------------
        self.tomo_time = 'To Estimate tomography time:'
        self.tomo_est_time = tk.Label(self.setupfrm, text=self.tomo_time)
        self.tomo_est_time.grid(row=row, column=self.col, columnspan=15, sticky="w")
        tomo_cal_button = tk.Button(self.setupfrm, text='Cal tomo time', command=self.cal_tomo_time)
        tomo_cal_button.grid(row=row,column=self.col+2, columnspan=2,sticky='w')
        #-------------------------xyl: now below open folder is reading the user folder to show, may no need below
        #self.pardir = tk.StringVar()
        #self.pardir.set(self.pvComm.getDir())
        #scdir_txt = tk.Button(
        #   self.setupfrm, text="Current data save directory:", command=self.getUserDir
        #)
        #scdir_txt.grid(row=row, column=self.col, sticky="w")
        #pardir_out = tk.Label(self.setupfrm, textvariable=self.pardir)
        #pardir_out.grid(row=row, column=self.col + 1, sticky="w", columnspan=3)
        #--------------------------------------------------------------------------------------------------------
        row += 2
        tot_est_txt = tk.Label(self.setupfrm, text="Estimated total scan time:")
        tot_est_txt.grid(row=row, column=self.col, sticky="w")
        self.tot_time = tk.StringVar()
        self.tot_time.set("0")
        tot_est_val = tk.Label(self.setupfrm, textvariable=self.tot_time)
        tot_est_val.grid(row=row, column=self.col + 1, sticky="w")

        #------------------------add scan box for xanes----------------------------------------------
        self.energy_step = 0.5
        self.energy_width = 100
        self.dwell_xanes = 20    #pre-defined value
        self.c_eng = 10
        
        self.xanes_frame = tk.Frame(self.setupfrm, borderwidth=1, relief="solid")
        self.xanes_frame.grid(row=row+3, column=self.col, columnspan=12, pady=2, padx=2, sticky="w")  #a sperate frame for xanes collection

        #xanes_paras = ['Energy (keV):','Energy step (eV):', 'Energy width (eV):', 'Dwell (ms):']
        xanes_paras_val = [self.c_eng, self.energy_step, self.energy_width, self.dwell_xanes]
        
        xanes_title = tk.Label(self.xanes_frame, text='XANES set up:')
        xanes_title.grid(row=row+3, column=self.col, sticky='w',pady=1,padx=2)  
        for xi, xp in enumerate(self.inputs_labels[2]):
            label = tk.Label(self.xanes_frame, text=xp)
            entry_var = tk.StringVar(value=f'{xanes_paras_val[xi]}')
            entry_var.trace("w", lambda name, index, mode, entry_var=entry_var: self.calcTime())
            entry = tk.Entry(self.xanes_frame, width=8,textvariable=entry_var)    
            label.grid(row=row+3, column=self.col+(xi%4)*2+1, sticky='w',pady=2)                               
            entry.grid(row=row+3, column=self.col+(xi%4)*2+2, sticky='w', padx=2, pady=2)
            self.scanParms.update({xp: entry})     # add values to scanParms  
        #-------------------------add a button to export fine images-----------------------------------
        
        self.export_frame = tk.Frame(self.setupfrm, borderwidth=1, relief="solid")
        self.export_frame.grid(row=51, column=self.col, columnspan=10, pady=1, padx=1, sticky="w")
        
        self.enter_value = tk.StringVar()
        self.enter_value_elm = tk.StringVar()
        scan_export_text = tk.Label(self.export_frame, text='scan to export:',width=14)
        scan_export_text.grid(row=51,column=10, padx=1,sticky='w')
        scan_range = tk.Entry(
                    self.export_frame,
                    width=10,
                    textvariable=self.enter_value)
        scan_range.grid(row=51,column=self.col+1, padx=1,sticky='w')
        scan_elm_text = tk.Label(self.export_frame, text='element:')
        scan_elm_text.grid(row=51,column=self.col+2, padx=1)
        scan_elm = tk.Entry(
                    self.export_frame,
                    width=3,
                    textvariable=self.enter_value_elm)
        scan_elm.grid(row=51,column=self.col+3, padx=1)# a entry for enter the scan start and end
        
        export_button = tk.Button(self.export_frame, text='export', command=self.export_scan_images)
        export_button.grid(row=51,column=self.col+4, padx=1)
        #self.export_prog_txt = tk.StringVar()
        self.export_prog_txt = "Done:"
        self.export_prog = tk.Label(self.export_frame, text=self.export_prog_txt)
        self.export_prog.grid(row=52, column=10,columnspan=10,padx=1,sticky='w')  #a text 
        note = tk.Label(self.export_frame, text='Batch export: [start,end,step] = [1,10] or [1,10,2]; 1 or >2 scans')
        note.config(fg='green')
        note.grid(row=53, column=10,columnspan=10,padx=1,sticky='w')  #a text 

        #------------------------try to get user folder when open GUI-----------------------------------
        #self.getUserDir()
        #self.choose_folder()
        
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
