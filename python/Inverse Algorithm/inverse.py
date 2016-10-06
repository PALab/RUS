# import Kivy
import kivy
import math
import random
import sys
import subprocess
import rus_parser as parser
import rus_forward as forward
import itertools
import graph1
from kivy.utils import get_color_from_hex as rgb



from kivy.app import App
from kivy.uix.button import Button
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.dropdown import DropDown
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.spinner import Spinner


shapeSpinner = Spinner(
    # default value shown
            text='Choose the Shape',
    # available values
            values=('Rectangle', 'Ellipsoidal Cylinder', 'Spheroid'),
    # just for positioning in our example
            width=150
            )
			
scSpinner = Spinner(
            text='Choose the number',
            values=('2', '3','5','6','9'),
			width=150
            )

hexSpinner = Spinner(
            text='Choose the type',
            values=('VTI', 'HTI'),
            width=150
            )
			
layout = GridLayout(cols=3,rows=17)

continueButton = Button(text='Continue')
sp = Button(text='Continue')
snLabel = Label(text='Enter the value of stiffness coefficient :', size_hint_x=None, width=150,)
hexLabel = Label(text='Hex Type :', size_hint_x=None, width=150,)
emptyCol= Label(text='', size_hint_x=None, width=150,)
emptyCol2 = Label(text='', size_hint_x=None, width=150,)
emptyCol3 = Label(text='', size_hint_x=None, width=150,)
emptyCol4 = Label(text='', size_hint_x=None, width=150,)
emptyCol5 = Label(text='', size_hint_x=None, width=150,)
emptyCol6 = Label(text='', size_hint_x=None, width=150)
emptyCol7 = Label(text='', size_hint_x=None, width=150)
emptyC12 = Label(text='X', size_hint_x=None, width=150)
emptyC13 = Label(text='X', size_hint_x=None, width=150)
emptyC11 = Label(text='X', size_hint_x=None, width=150)
emptyC66 = Label(text='X', size_hint_x=None, width=150)
emptyC22 = Label(text='X', size_hint_x=None, width=150)
emptyC23 = Label(text='X', size_hint_x=None, width=150)
emptyC33 = Label(text='X', size_hint_x=None, width=150)
emptyC55 = Label(text='X', size_hint_x=None, width=150)
c11 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c11')
c12 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c12')
c13 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c13')
c22 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c22')
c23 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c23')
c33 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c33')
c44 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c44')
c55 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c55')
c66 = TextInput(text='', multiline=False, size_hint_x=None,width=150, hint_text ='c66')
shape = '0'
sc = '0'
hextype = '0'
outputArray = []


class RUS(App):
# layout
    def build(self):
        
        sp.bind(on_press=self.spClicked)
        layout.add_widget(Label(text='Order of polynomials used to estimate the eigenvectors :', size_hint_x=None, width=400,))
        self.vectors = TextInput(text='6', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.vectors)
        layout.add_widget(Label(text='', size_hint_x=None, width=150))

        layout.add_widget(Label(text='Shape', size_hint_x=None, width=400))
        layout.add_widget(shapeSpinner)
        layout.add_widget(Label(text='', size_hint_x=None, width=150))

        layout.add_widget(Label(text='Dimension 1  :', size_hint_x=None, width=400,))
        self.dimension1 = TextInput(text='3.67', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.dimension1)
        layout.add_widget(Label(text='cm', size_hint_x=None, width=150,))

        layout.add_widget(Label(text='Dimension 2 :', size_hint_x=None, width=400,))
        self.dimension2 = TextInput(text='3.67', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.dimension2)
        layout.add_widget(Label(text='cm', size_hint_x=None, width=150))

        layout.add_widget(Label(text='Dimension 3 :', size_hint_x=None, width=400))
        self.dimension3 = TextInput(text='3.67', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.dimension3)
        layout.add_widget(Label(text='cm', size_hint_x=None, width=150))

        layout.add_widget(Label(text='Density :', size_hint_x=None, width=400))
        self.density = TextInput(text='1.7', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.density)
        layout.add_widget(Label(text='grams/cm^3', size_hint_x=None, width=150))
		
        layout.add_widget(Label(text='Number of stiffness coefficient :', size_hint_x=None, width=400))
        layout.add_widget(scSpinner)
        layout.add_widget(Label(text='', size_hint_x=None, width=150))
		
        layout.add_widget(Label(text='Minimum frequency  :', size_hint_x=None, width=400,))
        self.minfreq = TextInput(text='0.02', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.minfreq)
        layout.add_widget(Label(text='(set >1 KHz )', size_hint_x=None, width=150))
		
        layout.add_widget(Label(text='Maximum frequency  :', size_hint_x=None, width=400,))
        self.maxfreq = TextInput(text='0.11', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.maxfreq)
        layout.add_widget(Label(text='(set >5 or 10KHz)', size_hint_x=None, width=150))
		
        layout.add_widget(Label(text='Iterations  :', size_hint_x=None, width=400,))
        self.iterations = TextInput(text='1', multiline=False,size_hint_x=None, width=150)
        layout.add_widget(self.iterations)
        layout.add_widget(Label(text='', size_hint_x=None, width=150))
		
        layout.add_widget(snLabel)
        layout.add_widget(c11)
        layout.add_widget(c12)
        layout.add_widget(emptyCol)
        layout.add_widget(c13)
        layout.add_widget(c22)
        layout.add_widget(emptyCol2) 
        layout.add_widget(c23)
        layout.add_widget(c33)
        layout.add_widget(emptyCol3) 
        layout.add_widget(c44)
        layout.add_widget(c55)
        layout.add_widget(emptyCol4) 
        layout.add_widget(c66) 
        layout.add_widget(emptyCol5)
        layout.add_widget(sp)
       
        return layout
		
		
    @staticmethod	
    def remove_col():
        layout.remove_widget(c11)
        layout.remove_widget(c12)
        layout.remove_widget(c13)
        layout.remove_widget(c22)
        layout.remove_widget(c23)
        layout.remove_widget(c33)
        layout.remove_widget(c44)
        layout.remove_widget(c55)
        layout.remove_widget(c66)
        layout.remove_widget(snLabel)
        layout.remove_widget(continueButton)
        layout.remove_widget(sp)
        layout.remove_widget(hexLabel)
        layout.remove_widget(hexSpinner)
        layout.remove_widget(emptyCol)
        layout.remove_widget(emptyCol2)
        layout.remove_widget(emptyCol3)
        layout.remove_widget(emptyCol4)
        layout.remove_widget(emptyCol5)
        layout.remove_widget(emptyCol6)
        layout.remove_widget(emptyCol7)
        layout.remove_widget(emptyC12)
        layout.remove_widget(emptyC13)
        layout.remove_widget(emptyC22)
        layout.remove_widget(emptyC23)
        layout.remove_widget(emptyC33)
        layout.remove_widget(emptyC11)
        layout.remove_widget(emptyC55)
        layout.remove_widget(emptyC66)
		
    def show_selected_value(shapeSpinner, text):
        
        global shape
        global sc
        global hextype
        if text == 'Rectangle':
            shape = '0'
        if text == 'Ellipsoidal Cylinder':
            shape = '1'
        if text == 'Spheroid':
            shape = '2'
            
        if text == '2':
            sc = '2'
        if text == '3':
            sc = '3'
        if text == '5':
            sc = '5'
        if text == '6':
            sc = '6'
        if text == '9':
            sc = '9'
			
        if text == 'VTI':
            hextype = '1'
        if text == 'HTI':
            hextype = '2'
       
        if text == '5' and shape == '2':
            RUS.remove_col()
            layout.add_widget(hexLabel)
            layout.add_widget(hexSpinner)
            layout.add_widget(emptyCol)
            layout.add_widget(snLabel)
            layout.add_widget(emptyC11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol6)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '2' and shape == '2':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(emptyC12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(emptyC23)
            layout.add_widget(emptyC33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(emptyC66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '3' and shape == '2':	
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(emptyC23)
            layout.add_widget(emptyC33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(emptyC66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '6' and shape == '2':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '9' and shape == '2':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(c13)
            layout.add_widget(c22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(c55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '5':
            RUS.remove_col()
            layout.add_widget(hexLabel)
            layout.add_widget(hexSpinner)
            layout.add_widget(emptyCol)
            layout.add_widget(snLabel)
            layout.add_widget(emptyC11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol6)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '2':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(emptyC12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(emptyC23)
            layout.add_widget(emptyC33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(emptyC66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '3':	
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(emptyC23)
            layout.add_widget(emptyC33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(emptyC66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '6':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(emptyC13)
            layout.add_widget(emptyC22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(emptyC55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
        elif text == '9':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(c13)
            layout.add_widget(c22)
            layout.add_widget(emptyCol2) 
            layout.add_widget(c23)
            layout.add_widget(c33)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c44)
            layout.add_widget(c55)
            layout.add_widget(emptyCol4) 
            layout.add_widget(c66) 
            layout.add_widget(emptyCol5)
            layout.add_widget(sp)
       
    shapeSpinner.bind(text=show_selected_value)
    scSpinner.bind(text=show_selected_value)
    hexSpinner.bind(text=show_selected_value)
	
    @staticmethod
    def graphOutput():
        global outputArray
        output = open("freq_data", "r")
        lines = output.readlines()
        
        freq = []		
        yaxis = [] 
        predicted = []
        predictedCount = 0
        splitFirst = 0
        splitSecond = 15
        
        for i in range(0,len(lines)):
            outputArray.append(lines[i])
                   
        for i in range(0,len(lines)):
            outputArray[i] = (outputArray[i])
            
            #1 indicated to split no more than 1 times, and [0] means take the first element returned by split
            splitData = outputArray[i].split("\n",1)[0]
            xValue = splitData.split("\t",1)[0]
            yValue = splitData.split("\t")[-1]
            yaxis.append(yValue)
            freq.append(xValue) 
            
        predictedFreq = open("predictedf", "r")
        predictedLines = predictedFreq.readlines()
        
        
        while (predictedCount < int(yaxis[0])):
            firstP = ''.join(predictedLines)[splitFirst:splitSecond]
            if firstP[-1:] == '.' :
                firstP = ''.join(predictedLines)[splitFirst:splitSecond-1 ] # predictedFreq file doesn't always have 13 d.p (sometimes only 12 d.p)
            predicted.append(firstP)
            predictedCount = predictedCount + 1
            splitFirst = splitFirst + 15
            splitSecond = splitSecond + 15
        
        colors = itertools.cycle([
            rgb('7dac9f'), rgb('dc7062'), rgb('66a8d4'), rgb('e5b060')])
        graph_theme = {
            'label_options': {
                'color': rgb('444444'),  # color of tick labels and titles
                'bold': True},
            'background_color': rgb('f8f8f2'),  # back ground color of canvas
            'tick_color': rgb('808080'),  # ticks and grid
            'border_color': rgb('808080')}  # border drawn around each graph

        graph = graph1.Graph(
            xlabel='Frequency (Hz) x0.001',
            ylabel='Unit',
            x_ticks_major=5,
            y_ticks_major=1,
            y_grid_label=True,
            x_grid_label=True,
            padding=5,
            xlog=False,
            ylog=False,
            x_grid=False,
            y_grid=True,
            #assume that freq_data start in the lowest freq
            xmin=int(math.fsum([float(freq[1])] * 1000)) -5,
            xmax=int(math.fsum([float(freq[len(freq)-14])] * 1000)),
            ymin=0,
            ymax=1,
            draw_border = True,
            **graph_theme)
				
        count = 0
        y = 1
        
        # Display the number of freq_data (stated in line 1)
        while(count < int(yaxis[0])):
            if(count == 0):
                plot = graph1.SmoothLinePlot(color=rgb('7dac9f'))
                plot.points = [(count, x) for x in range(0, 0)]
                graph.add_plot(plot)
                count = count + 1
            else:		
                freqx100 = math.fsum([float(freq[y])] * 1000)
                yaxisToInt = math.fsum([float(yaxis[y])] )
                plot = graph1.SmoothLinePlot(color=rgb('7dac9f'))
                plot.points = [(int(freqx100), x) for x in range(0,int(yaxisToInt)+1)]
                graph.add_plot(plot)
                count = count + 1
                y = y + 1
            
        
        # Display the freq_data file
        while(count < len(freq)-14):
            if(count == 0):
                y = 0
                plot = graph1.SmoothLinePlot(color=rgb('dc7062'))
                plot.points = [(count, x) for x in range(0, 0)]
                graph.add_plot(plot)
                count = count + 1
            else:		
                freqx100 = math.fsum([float(freq[y])] * 1000)
                yaxisToInt = math.fsum([float(yaxis[y])] )
                plot = graph1.SmoothLinePlot(color=rgb('dc7062'))
                plot.points = [(int(freqx100), x) for x in range(0,int(yaxisToInt)+1)]
                graph.add_plot(plot)
                count = count + 1
                y = y + 1
        
        count = 0
        y = 0
        # Display the predictedFreq file
        while(count < int(yaxis[0])):
            if(count == 0):
                plot = graph1.SmoothLinePlot(color=rgb('e5b060'))
                plot.points = [(count, x) for x in range(0, 0)]
                graph.add_plot(plot)
                count = count + 1
            else:		
                freqx100 = math.fsum([float(predicted[y])] * 1000)
                plot = graph1.SmoothLinePlot(color=rgb('e5b060'))
                plot.points = [(int(freqx100), x) for x in range(0,2)]
                graph.add_plot(plot)
                count = count + 1
                y = y + 1
        
        
        layout.add_widget(graph)
		
    

# button click function
    def spClicked(self,btn):
	
        if shape == '2':
			if sc == '2':	
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' + ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text + ' --c44 ' + c44.text
			elif sc == '3':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' +  ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c12 ' + c12.text + ' --c44 ' + c44.text
			elif sc == '5' and hextype == '1':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' + ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c33 ' + c33.text +  ' --c23 ' + c23.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text + ' --c66 ' + c66.text + ' --hextype ' + hextype
			elif sc == '5' and hextype == '2':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' + ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c33.text +  ' --c33 ' + c23.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text + ' --c66 ' + c66.text + ' --hextype ' + hextype
			elif sc == '6':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' +  ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c33 ' + c33.text + ' --c23 ' + c23.text + ' --c12 ' + c12.text +  ' --c44 ' + c44.text + ' --c66 ' + c66.text
			elif sc == '9':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 1' +  ' --d3 1' +  ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c22 ' + c22.text + ' --c33 ' + c33.text + ' --c23 ' + c23.text +  ' --c13 ' + c13.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text +  ' --c55 ' + c55.text + ' --c66 ' + c66.text
        else:
			if sc == '2':	
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text + ' --c44 ' + c44.text
			elif sc == '3':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c12 ' + c12.text + ' --c44 ' + c44.text
			elif sc == '5' and hextype == '1':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c33 ' + c33.text +  ' --c23 ' + c23.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text + ' --c66 ' + c66.text + ' --hextype ' + hextype
			elif sc == '5' and hextype == '2':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c33.text +  ' --c33 ' + c23.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text + ' --c66 ' + c66.text + ' --hextype ' + hextype
			elif sc == '6':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c33 ' + c33.text + ' --c23 ' + c23.text + ' --c12 ' + c12.text +  ' --c44 ' + c44.text + ' --c66 ' + c66.text
			elif sc == '9':
				command = 'python rus.py inverse --order ' + self.vectors.text + ' --freqmin ' + self.minfreq.text + ' --freqmax ' + self.maxfreq.text + ' --iteration ' + self.iterations.text + ' --shape ' + shape + ' --d1 ' + self.dimension1.text + ' --d2 ' + self.dimension2.text + ' --d3 ' + self.dimension3.text + ' --rho ' + self.density.text + ' --ns ' + sc + ' --c11 ' + c11.text +  ' --c22 ' + c22.text + ' --c33 ' + c33.text + ' --c23 ' + c23.text +  ' --c13 ' + c13.text + ' --c12 ' + c12.text + ' --c44 ' + c44.text +  ' --c55 ' + c55.text + ' --c66 ' + c66.text
        
        print (command)
        layout.clear_widgets()
        f = open("output.txt", "w")
        subprocess.call(command, shell=True , stdout=f)
        output = open("output.txt", "r")
        readOutput = output.readlines()
        
        if sc == '2':
            label = Label(text = 'Red lines = Freq_Data\n' + 'Blue lines = Starting\n' + 'Orange lines = Ending\n\n'+
                                readOutput[len(readOutput)-8] + ('C11 = ' + readOutput[len(readOutput)-5]) + ('C44 = ' + readOutput[len(readOutput)-4]), size_hint_x=None, width=200)
        if sc == '3':
            label = Label(text = 'Red lines = Freq_Data\n' + 'Blue lines = Starting\n' + 'Orange lines = Ending\n\n'+
                                readOutput[len(readOutput)-9] + ('C11 = ' + readOutput[len(readOutput)-6]) + ('C12 = ' + readOutput[len(readOutput)-5])+ ('C44 = ' + readOutput[len(readOutput)-4]),size_hint_x=None, width=200)
        if sc == '5':
            label = Label(text = 'Red lines = Freq_Data\n' + 'Blue lines = Starting\n' + 'Orange lines = Ending\n\n'+
                                readOutput[len(readOutput)-11] + ('C12 = ' + readOutput[len(readOutput)-8]) + ('C23 = ' + readOutput[len(readOutput)-7])+ ('C33 = ' + readOutput[len(readOutput)-6]) + ('C44 = ' + readOutput[len(readOutput)-5]) + ('C66 = ' + readOutput[len(readOutput)-4]),size_hint_x=None, width=200)
        if sc == '6':
            label = Label(text = 'Red lines = Freq_Data\n' + 'Blue lines = Starting\n' + 'Orange lines = Ending\n\n'+
                                readOutput[len(readOutput)-12] + ('C11 = ' + readOutput[len(readOutput)-9]) + ('C12 = ' + readOutput[len(readOutput)-8])+ ('C23 = ' + readOutput[len(readOutput)-7]) + ('C33 = ' + readOutput[len(readOutput)-6]) + ('C44 = ' + readOutput[len(readOutput)-5]) + ('C66 = ' + readOutput[len(readOutput)-4]) ,size_hint_x=None, width=200)
        if sc == '9':
            label = Label(text = 'Red lines = Freq_Data\n' + 'Blue lines = Starting\n' + 'Orange lines = Ending\n\n'+
                                readOutput[len(readOutput)-15] + ('C11 = ' + readOutput[len(readOutput)-12]) + ('C12 = ' + readOutput[len(readOutput)-11])+ ('C13 = ' + readOutput[len(readOutput)-10]) + ('C22 = ' + readOutput[len(readOutput)-9]) + ('C23 = ' + readOutput[len(readOutput)-8]) + ('C33 = ' + readOutput[len(readOutput)-7]) + ('C44 = ' + readOutput[len(readOutput)-6]) + ('C55 = ' + readOutput[len(readOutput)-5]) + ('C66 = ' + readOutput[len(readOutput)-4]),size_hint_x=None, width=200)
        layout.add_widget(label)  
        RUS.graphOutput()
        
# run app
if __name__ == "__main__":
    RUS().run()

 