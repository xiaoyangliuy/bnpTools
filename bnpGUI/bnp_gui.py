#!/home/beams/USERBNP/.conda/envs/py36/bin/python

import tkinter as tk
from tkinter import ttk
from setupFrame import setupFrame
from scanFrame import scanFrame
from tkinter import messagebox
from epics import caget, caput
import os
import shutil
from tkinter import font
#-------------------add a window for input user and scan id-----------------------
class UserinfoWindow(tk.Toplevel):
    def __init__(self,parent):
        super().__init__(parent)
        self.title('Input')
        self.geometry('800x600')
        #input_window = tk.Toplevel()
        self.inp_frm = tk.Frame(self, borderwidth=1, relief="solid", width=600, height=400,padx=10,pady=10)
        self.inp_frm.grid(row=3, column=3, sticky="ew",columnspan=3)
        font_size_1 = font.Font(size=16)
        font_size_2 = font.Font(size=14)
        #--------------------PVs' value----------------------------------------------
        self.cyclefolder = caget('9idbBNP:saveData_fileSystem')  #root till cycle
        #self.cyclefolder = self.cyclefolder.replace('//micdata/data1','/mnt/micdata1')
        self.mda_dir = caget('9idbBNP:saveData_subDir')
        self.scan_num = caget('9idbBNP:saveData_scanNumber')
        self.det_file = caget('9idbXMAP:netCDF1:FilePath', as_string=True)
        #--------------------------Pv values read from current PV-------------------------------------------------
        self.cycle_label = tk.Label(self.inp_frm,text='Cycle:',width=12,font=font_size_1)
        self.cycle_label.grid(row=1, column=1, sticky="w")
        
        self.cycle_text = tk.StringVar(value=self.cyclefolder)
        self.cycle_entry = tk.Entry(self.inp_frm,textvariable=self.cycle_text,width=38,font=font_size_2)
        self.cycle_entry.grid(row=1, column=2, sticky="w")
        
        self.user_label = tk.Label(self.inp_frm,text="User:",width=12, font=font_size_1)
        self.user_label.grid(row=2, column=1, sticky="w")

        self.user_text = tk.StringVar(value=self.mda_dir)
        self.user_entry = tk.Entry(self.inp_frm,textvariable=self.user_text,width=38,font=font_size_2)
        self.user_entry.grid(row=2, column=2, sticky="w")
        
        self.scan_label = tk.Label(self.inp_frm, text="Scan start:",width=12,font=font_size_1)
        self.scan_label.grid(row=3, column=1, sticky="w")
        
        self.scan_text = tk.StringVar(value=self.scan_num)
        self.scan_entry = tk.Entry(self.inp_frm,textvariable=self.scan_text,width=38,font=font_size_2)
        self.scan_entry.grid(row=3, column=2, sticky="w")
        
        self.det_label = tk.Label(self.inp_frm,text="Det:", width=12,font=font_size_1)
        self.det_label.grid(row=4, column=1, sticky="w")
        
        self.det_text = tk.StringVar(value=self.det_file)
        self.det_entry = tk.Entry(self.inp_frm,textvariable=self.det_text,width=38,font=font_size_2)
        self.det_entry.grid(row=4, column=2, sticky="w")
        #---------------------------------------------------------------------------
        self.submit_button = tk.Button(self.inp_frm, text="Create", command=self.submit_input,width=12,font=font_size_1)
        self.submit_button.grid(row=5, column=1, sticky="ew",padx=2,pady=2)
        
        self.ok_button = tk.Button(self.inp_frm,text='OK',command=self.accept,width=12,font=font_size_1)
        self.ok_button.grid(row=5,column=2,sticky='ew',padx=2,pady=2)
    
    def submit_input(self):
        cycle = self.cycle_entry.get() #BNP+cycle
        user = self.user_entry.get()  #username+mda
        scan = self.scan_entry.get()    #next scan ID        
        det = self.det_entry.get()   #username+flyXRF
        user_name = user.split('/')[0]+'/'
        #-------------------------make mda and flyXRF file in the user folder-----------
        cycle2 = cycle.replace('//micdata/data1','/mnt/micdata1')
        os.makedirs(f'{cycle2}/{user}')
        user_det = user.replace('mda','flyXRF')
        os.makedirs(f'{cycle2}/{user_det}')
        #--------------------------------------------------------------------------------
        #copy override and standard files from a fixed folder in current cycle----------
        fixfolder = f'{cycle2}/test_gui_xyl3'
        override = 'maps_fit_parameters_override.txt'
        std = 'maps_standardinfo.txt'       
        source_override = os.path.join(fixfolder,override)
        source_std = os.path.join(fixfolder,std)
        des_override = os.path.join(cycle2,user_name)
        des_std = os.path.join(cycle2,user_name)       
        shutil.copy2(source_override, des_override)
        shutil.copy2(source_std, des_std)
        #---------------------------change PVs values------------------------------
        caput('9idbBNP:saveData_fileSystem',cycle,wait=True)
        caput('9idbBNP:saveData_subDir',user,wait=True)
        caput('9idbBNP:saveData_scanNumber',scan,wait=True)
        caput('9idbXMAP:netCDF1:FilePath',det,wait=True)
        
        self.destroy()
        
    def accept(self):
        self.destroy()
        #self.parent.deiconify()
        
        

class MainWindow(tk.Frame):
        
    def __init__(self, master):        
        
        self.master = master
        master.title("BNP Scan")
        self.show_user_info()
        #self.user_info_window = UserinfoWindow(self.master)
        #self.wait_window(self.user_info_window)


        self.tabControl = ttk.Notebook(master)
        setup_tab = setupFrame(self.tabControl)
        scan_tab = scanFrame(self.tabControl, setup_tab)
        setup_tab.addToScanBtn(scan_tab.insertScan)

        self.tabControl.add(setup_tab.setupfrm, text ='Scan Setup')
        self.tabControl.add(scan_tab.scanfrm, text ='Data Collection')
        self.tabControl.pack(expand=1, fill="both")
        
        #self.show_user_info()
        

    def show_user_info(self):
        user_info_window = UserinfoWindow(self.master)
        self.master.wait_window(user_info_window)    
if __name__ == "__main__":
    #-----------------------add window for input user info and scan id-----------------
    #input_window = tk.Toplevel()
    #inp = UserinfoWindow(input_window)
   # inp.geometry('860x680')
    #input_window.mainloop()
    
    #-------------------------------------------------------------------------
    root = tk.Tk()
    

    #user_info_window = UserinfoWindow(root)
    #root.wait_window(user_info_window)
    gui = MainWindow(root)
    #gui.show_user_info()
    root.geometry('1550x780')
    root.resizable(True, True)
    root.mainloop()

