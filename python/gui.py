# import Kivy
import kivy
import random

from kivy.app import App
from kivy.uix.button import Button
from kivy.uix.gridlayout import GridLayout
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
            width=80
            )
			
scSpinner = Spinner(
            text='Choose the number',
            values=('2', '3','5','6','9'),
			width=80
            )

hexSpinner = Spinner(
            text='Choose the type',
            values=('1.VTI', '2.HTI'),
            width=80
            )
			
layout = GridLayout(cols=3,rows=14,row_default_height=50)

continueButton = Button(text='Continue')
snLabel = Label(text='Enter the value of stiffness coefficient :', size_hint_x=None, width=500,)
hexLabel = Label(text='Hex Type :', size_hint_x=None, width=500,)
emptyCol= Label(text='', size_hint_x=None, width=80,)
emptyCol2 = Label(text='', size_hint_x=None, width=80,)
emptyCol3 = Label(text='', size_hint_x=None, width=80,)
emptyCol4 = Label(text='', size_hint_x=None, width=80,)
emptyCol5 = Label(text='', size_hint_x=None, width=80,)
c11 = TextInput(text='4.284', multiline=False, size_hint_x=None,width=80)
c12 = TextInput(text='4.284', multiline=False, size_hint_x=None,width=80)
c13 = TextInput(text='4.284', multiline=False, size_hint_x=None,width=80)
c22 = TextInput(text='2.448', multiline=False, size_hint_x=None,width=80)
c23 = TextInput(text='4', multiline=False, size_hint_x=None,width=80)
c33 = TextInput(text='2.448', multiline=False, size_hint_x=None,width=80)
c44 = TextInput(text='2.448', multiline=False, size_hint_x=None,width=80)
c55 = TextInput(text='5.508', multiline=False, size_hint_x=None,width=80)
c66 = TextInput(text='1.7', multiline=False, size_hint_x=None,width=80)
class RUS(App):
# layout
    def build(self):
        
        '''	
#  drop down box for shape	  
        
        shapeDrop = DropDown()
        shapes = ['1.Rectangle', '2.Ellipsoidal Cylinder', '3.Spheroid']
        for shape in shapes:
            btn = Button(text='%r' % shape, size_hint_y=None, height=30)
            btn.bind(on_release=lambda btn: shapeDrop.select(btn.text))
            shapeDrop.add_widget(btn)
        mainbutton = Button(text='Choose the shape', size_hint=(1, 1))
        mainbutton.bind(on_release=shapeDrop.open)
        shapeDrop.bind(on_select=lambda instance, x: setattr(mainbutton, 'text', x))
#  drop down box for hex type		
        hexDrop = DropDown()
        hex = ['1.VTI', '2.HTI']
        for type in hex:
            btn1 = Button(text='%r' % type, size_hint_y=None, height=30)
            btn1.bind(on_release=lambda btn1: hexDrop.select(btn1.text))
            hexDrop.add_widget(btn1)
        hexbutton = Button(text='Hex Type', size_hint=(1, 1))
        hexbutton.bind(on_release=hexDrop.open)
        hexDrop.bind(on_select=lambda instance, x: setattr(hexbutton, 'text', x))
#  drop down box for stiffness coefficient		
        stiffnessDrop = DropDown()
        sitffness = ['2', '3','5','6','9']
        for i in sitffness:
            btn2 = Button(text='%r' % i, size_hint_y=None, height=30)
            btn2.bind(on_release=lambda btn2: stiffnessDrop.select(btn2.text))
            stiffnessDrop.add_widget(btn2)
        sfbutton = Button(text='Stiffness', size_hint=(1, 1))
        sfbutton.bind(on_release=stiffnessDrop.open)
        stiffnessDrop.bind(on_select=lambda instance, x: setattr(sfbutton, 'text', x))
        '''
        layout.add_widget(Label(text='Number of eigen frequencies to print :', size_hint_x=None, width=500))
        layout.add_widget(TextInput(text='10', multiline=False,size_hint_x=None, width=150))
        layout.add_widget(Label(text='', size_hint_x=None, width=80))

        layout.add_widget(Label(text='Order of the polynomial to use to estimate eigenvectors :', size_hint_x=None, width=100,))
        layout.add_widget(TextInput(text='6', multiline=False,size_hint_x=None, width=80))
        layout.add_widget(Label(text='', size_hint_x=None, width=100,))

        layout.add_widget(Label(text='Shape', size_hint_x=None, width=100,))
        layout.add_widget(shapeSpinner)
        layout.add_widget(Label(text='', size_hint_x=None, width=100,))

        layout.add_widget(Label(text='Dimension 1  :', size_hint_x=None, width=100,))
        layout.add_widget(TextInput(text='3.67', multiline=False,size_hint_x=None, width=80))
        layout.add_widget(Label(text='cm', size_hint_x=None, width=100,))

        layout.add_widget(Label(text='Dimension 2 :', size_hint_x=None, width=100,))
        layout.add_widget(TextInput(text='3.67', multiline=False,size_hint_x=None, width=80))
        layout.add_widget(Label(text='cm', size_hint_x=None, width=100,))

        layout.add_widget(Label(text='Dimension 3 :', size_hint_x=None, width=100,))
        layout.add_widget(TextInput(text='3.67', multiline=False,size_hint_x=None, width=80))
        layout.add_widget(Label(text='cm', size_hint_x=None, width=100,))

        layout.add_widget(Label(text='Density :', size_hint_x=None, width=100,))
        layout.add_widget(TextInput(text='1.7', multiline=False,size_hint_x=None, width=80))
        layout.add_widget(Label(text='grams/cm^3', size_hint_x=None, width=100,))
		
        layout.add_widget(Label(text='Number of stiffness coefficient :', size_hint_x=None, width=100,))
        layout.add_widget(scSpinner)
        layout.add_widget(Label(text='', size_hint_x=None, width=100,))
        layout.add_widget(continueButton)
       



        
        
        return layout
		
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
        layout.remove_widget(hexLabel)
        layout.remove_widget(hexSpinner)
        layout.remove_widget(emptyCol)
        layout.remove_widget(emptyCol2)
        layout.remove_widget(emptyCol3)
        layout.remove_widget(emptyCol4)
        layout.remove_widget(emptyCol5)
		
    def show_selected_value(scSpinner, text):
        
        if text == '5':
            RUS.remove_col()
            layout.add_widget(hexLabel)
            layout.add_widget(hexSpinner)
            layout.add_widget(emptyCol)
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol2)
            layout.add_widget(c13)
            layout.add_widget(c22)
            layout.add_widget(emptyCol3) 
            layout.add_widget(c23)
            layout.add_widget(emptyCol4)
            layout.add_widget(continueButton)
        elif text == '2':
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(continueButton)
            
        elif text == '3':	
            RUS.remove_col()
            layout.add_widget(snLabel)
            layout.add_widget(c11)
            layout.add_widget(c12)
            layout.add_widget(emptyCol)
            layout.add_widget(c13)
            layout.add_widget(emptyCol2)
            layout.add_widget(continueButton)
        elif text == '6':
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
            layout.add_widget(continueButton)
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
            layout.add_widget(continueButton)
        
    scSpinner.bind(text=show_selected_value)
		
    

# button click function
    #def buttonClicked(self,btn):
        #self.lbl1.text = "You wrote " + self.txt1.text

# run app
if __name__ == "__main__":
    RUS().run()
 