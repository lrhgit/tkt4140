"""
Example on generic model.py file which inspects the arguments
of the compute function and automatically generates a relevant
InputForm class.
"""

import wtforms

from compute import compute_Manufactured_solution as compute
import inspect

arg_names = inspect.getargspec(compute).args
defaults  = inspect.getargspec(compute).defaults
print "Follow the link and add 'MMS'"
class InputForm(wtforms.Form):
    pass

# Augment defaults with None elements for the positional
# arguments
defaults = [None]*(len(arg_names)-len(defaults)) + list(defaults)
# Map type of default to right form field
type2form = {type(1.0): wtforms.FloatField,
             type(1):   wtforms.IntegerField,
             type(''):  wtforms.TextField,
             }

for name, value in zip(arg_names, defaults):
    #print name, value
    if value is None:
        setattr(InputForm, name, wtforms.TextField(
            validators=[wtforms.validators.InputRequired()]))
    else:
        if type(value) in type2form:
            setattr(InputForm, name, type2form[type(value)](
                default=value,
                validators=[wtforms.validators.InputRequired()]))
        else:
            raise TypeError('argument %s %s not supported' %
                            name, type(value))

if __name__ == '__main__':
    for item in dir(InputForm):
        if item in arg_names:
            print item, getattr(InputForm, item)
