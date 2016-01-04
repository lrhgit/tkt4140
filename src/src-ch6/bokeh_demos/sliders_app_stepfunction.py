import logging
logging.basicConfig(level=logging.DEBUG)

import numpy as np
from bokeh.plotting import figure
from bokeh.models import Plot, ColumnDataSource
from bokeh.properties import Instance
from bokeh.server.app import bokeh_app
from bokeh.server.utils.plugins import object_page
from bokeh.models.widgets import HBox, Slider, TextInput, Select, VBoxForm


class SlidersApp(HBox):
    """An example of a browser-based, interactive plot with slider controls."""

    extra_generated_classes = [["SlidersApp", "SlidersApp", "HBox"]]
    inputs = Instance(VBoxForm)

    text = Instance(Select)
    solver = Instance(Select)
    CFL = Instance(Slider)
    time = Instance(Slider)
    Period = Instance(Slider)
    plot = Instance(Plot)
    source = Instance(ColumnDataSource)

    @classmethod
    def create(cls):
        """One-time creation of app's objects.

        This function is called once, and is responsible for
        creating all objects (plots, datasources, etc)
        """
        obj = cls()
        obj.source = ColumnDataSource(data=dict(x=[], y=[]))

        obj.text = Select(
            title="Select initial shape", name='text', value='step', options=["step", "sine"]
        )
        obj.solver = Select(
            title="Select solver", name='solver', value='ftbs', options=["ftbs", "macCormack", "lax_friedrich"]
        )

        obj.CFL = Slider(
            title="CFL", name='CFL',
            value=0.4, start=0.1, end=1.0, step=0.1
        )
        obj.time = Slider(
            title="time", name='time',
            value=0.0, start=0.0, end=1.0, step=0.1
        )
        obj.Period = Slider(
            title="Period T (applicable for sine)", name='Period',
            value=0.6, start=0.2, end=1.0, step=0.2
        )


        toolset = "crosshair,pan,reset,resize,save,wheel_zoom"
        # Generate a figure container
        plot = figure(title_text_font_size="12pt",
                      plot_height=400,
                      plot_width=400,
                      tools=toolset,
                      title=obj.text.value,
                      x_range=[0, 1.5],
                      y_range=[-0.25, 1.5]
        )

        # Plot the line by the x,y values in the source property
        plot.line('x', 'y', source=obj.source,
                  line_width=3,
                  line_alpha=0.6
        )

        obj.plot = plot
        obj.update_data()

        obj.inputs = VBoxForm(
            children=[
                obj.text, obj.solver, obj.CFL, obj.time, obj.Period
            ]
        )

        obj.children.append(obj.inputs)
        obj.children.append(obj.plot)

        return obj

    def setup_events(self):
        """Attaches the on_change event to the value property of the widget.

        The callback is set to the input_change method of this app.
        """
        super(SlidersApp, self).setup_events()
        if not self.text:
            return

        # Text box event registration
        self.text.on_change('value', self, 'input_change')

        # Slider event registration
        for w in ["solver", "CFL", "time", "Period"]:
            getattr(self, w).on_change('value', self, 'input_change')

    def input_change(self, obj, attrname, old, new):
        """Executes whenever the input form changes.

        It is responsible for updating the plot, or anything else you want.

        Args:
            obj : the object that changed
            attrname : the attr that changed
            old : old value of attr
            new : new value of attr
        """
        self.update_data()
        self.plot.title = self.text.value

    def update_data(self):
        """Called each time that any watched property changes.

        This updates the womersleyprofile, with the most recent change. the u value is normalized to one. I. e. the magnitude
        of the pressure gradient is not incorporated in this plot.
        This is stored as two numpy arrays in a dict into the app's
        data source property.
        """

        # Get the current slider values
        
        func = self.text.value
        CFL = self.CFL.value
        solver = self.solver.value
        time = self.time.value
        CFL = str(float(CFL))
        Period = self.Period.value
        Period = str(float(Period))
        time = str(float(time))
        print "CFL, time: ", CFL, time
        print solver
        
        if func =='step':
            X = np.array([])
            Y = np.array([])
            cpath = 'solutions/' + func +'/' + solver +'/' + CFL + '/'
            #with open(path + '0')
            tpath = cpath + time + '.txt'
            with open(tpath,'r') as filename:
                for line in filename:
                    Xtemp = float(line.split()[0])
                    Ytemp = float(line.split()[1])
                    X = np.append(X, Xtemp)
                    Y = np.append(Y, Ytemp)
            filename.close()
        elif func=='sine':
            X = np.array([])
            Y = np.array([])
            cpath = 'solutions/' + func +'/' + solver +'/' + Period +'/' + CFL + '/'
            #with open(path + '0')
            tpath = cpath + time + '.txt'
            with open(tpath,'r') as filename:
                for line in filename:
                    Xtemp = float(line.split()[0])
                    Ytemp = float(line.split()[1])
                    X = np.append(X, Xtemp)
                    Y = np.append(Y, Ytemp)
            filename.close()
            

        logging.debug(
            "PARAMS: CFL: %s time: %s", self.CFL.value,
            self.time.value
        )

        self.source.data = dict(x=X, y=Y)


# The following code adds a "/bokeh/sliders/" url to the bokeh-server. This
# URL will render this sliders app. If you don't want to serve this
# applet from a Bokeh server (for instance if you are embedding in a separate
# Flask application), then just remove this block of code.
@bokeh_app.route("/bokeh/sliders/")
@object_page("sin")
def make_sliders():
    app = SlidersApp.create()
    return app

