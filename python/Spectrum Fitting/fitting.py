import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import csv

import sys
import peakutils
from peakutils.plot import plot as pplot

from peak_picker import Peak, PeakPicker

if sys.version_info[0] < 3:
    from Tkinter import * #python 2

	
else:
    from tkinter import * #python 3
import tkMessageBox

class toolsbar(NavigationToolbar2TkAgg):
	def getmessage(self):
		return self.message
class fitting:
	
	f = None
	graphFrame = None
	canvas = None
	rightFrame = None
	xDate = None
	yDate = None
	peaks = None
	listbox = None
	mode = False
	picker = None
	press = None
	
	def modeSelect(self):
		self.mode = True
		self.toolbar.pack_forget()
		
		
	def modeSelect2(self):
		self.mode = False
		self.toolbar.pack(side=LEFT)
		
		
	def saveToFile(self):
		p = self.listbox.get(0, END)
		f = open("freq_data","w")
		for i in p:
			f.write("%s1\n" % str(i).ljust(24))

	

	def __init__(self, master):

		master.minsize(width=800, height=600)

		leftFrame = Frame(master, background="red")
		self.rightFrame = Frame(master, background="orange")
		rightFrame = self.rightFrame
		
		lefttopFrame = Frame(leftFrame, background="blue")
		lefttopFrame.pack(side=TOP)

		leftbotFrame = Frame(leftFrame, background="blue")
		leftbotFrame.pack(side=BOTTOM, fill=BOTH)

		scrollbar = Scrollbar(lefttopFrame, orient=VERTICAL)
		self.listbox = Listbox(lefttopFrame, yscrollcommand=scrollbar.set)
		self.listbox.config(height = 30)

		scrollbar.config(command=self.listbox.yview)
		scrollbar.pack(side=RIGHT, fill=BOTH)
		self.listbox.pack(side=TOP, fill=BOTH)

		

		

		leftFrame.pack(side = LEFT, fill=BOTH)
		rightFrame.pack(side = RIGHT, fill=BOTH)


		self.graphFrame = Frame(rightFrame, background="yellow")

		self.graphFrame.pack(pady = 0, padx = 0)
		
		self.save_peaks = Button(leftFrame, command=self.saveToFile, text="save peaks")
		self.save_peaks.pack(fill=BOTH)
		label = Label
		buttonsFrame = Frame(leftbotFrame, background="green")
		buttonsFrame.pack(side=BOTTOM)
		

	
		




		picker = PeakPicker(X_DATA, Y_DATA, Trigger_Point)
		self.picker = picker
		
	
		self.canvas = FigureCanvasTkAgg(picker.fig, root)
		self.canvas.show()

		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

		
		
		
		self.toolbar = toolsbar(self.canvas, root)

		self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

		self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
		

		self.modeSelect2 = Button(buttonsFrame, width=12, command=self.modeSelect2, text="Toolbar")
		self.modeSelect2.pack(side=BOTTOM, fill=Y)

		self.modeSelect = Button(buttonsFrame, width=12, command=self.modeSelect, text="Select Peaks")
		self.modeSelect.pack(side=BOTTOM, fill=Y)

		self.modelabel = Label(buttonsFrame, text = "Mode select:")
		self.modelabel.pack(side=BOTTOM, fill=Y)
		
		
		def canvas_mouse_click(event):
			if self.mode is True:
				picker.onclick(event)
				self.listbox.delete(0, END)
				for i in picker.peaks:	
					self.listbox.insert(END, i.center)
		self.canvas.mpl_connect('button_press_event' , canvas_mouse_click)
		
	


		


def on_closing():
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
	root.quit()
        root.destroy()


		
if __name__ == '__main__':
    with open("80-90.CSV") as csvfile:
    	readCSV = csv.reader(csvfile, delimiter=',')
	X_DATA = []
	Y_DATA = []
	Trigger_Point = 0
    	for row in readCSV:
		X_DATA.append(float(row[3]))
		Y_DATA.append(float(row[4]))
		if row[0] == "Trigger Point":
			Trigger_Point = float(row[1])
    X_DATA = np.asarray(X_DATA)
    Y_DATA = np.asarray(Y_DATA)


root = Tk()
a = fitting(root)

root.protocol("WM_DELETE_WINDOW", on_closing)

root.mainloop()
root.quit() 

