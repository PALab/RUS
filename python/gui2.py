import itertools
import graph1
import rus_parser as parser
import rus_forward as forward
from math import sin, cos, pi
from kivy.utils import get_color_from_hex as rgb
from kivy.uix.boxlayout import BoxLayout
from kivy.app import App
	
args = parser.start()
return_list=[]
if args.subcommand == 'forward':
    return_list=forward.start(args)
    
    
class RUSApp(App):

        def build(self):
            b = BoxLayout(orientation='vertical')
            # example of a custom theme
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
                xlabel='Points',
                ylabel='Frequency (kHz)',
                x_ticks_major=1,
                y_ticks_major=25,
                y_grid_label=True,
                x_grid_label=False,
                padding=5,
                xlog=False,
                ylog=False,
                x_grid=False,
                y_grid=True,
                xmin=0,
                xmax=len(return_list)+1,
                ymin=0,
                ymax=return_list[len(return_list)-1]/100,
                **graph_theme)
				
            
            count = 0
            y = 0
            print(len(return_list))
            
            while(count < len(return_list)+1):
                if(count == 0):
                    plot = graph1.SmoothLinePlot(color=next(colors))
                    plot.points = [(count, x) for x in range(0, 0)]
             # for efficiency, the x range matches xmin, xmax 
                    graph.add_plot(plot)
                    count = count + 1
                else:				
                    plot = graph1.SmoothLinePlot(color=next(colors))
                    plot.points = [(count, x) for x in range(0, int(return_list[y]/100))]
             # for efficiency, the x range matches xmin, xmax 
                    graph.add_plot(plot)
                    count = count + 1
                    y = y + 1
                
            b.add_widget(graph)
            return b
if __name__ == '__main__':
    RUSApp().run()