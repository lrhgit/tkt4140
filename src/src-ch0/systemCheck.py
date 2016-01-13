# src-ch0/systemCheck.py

def systemCheck():
    '''
    Check for necessary moduls needed in the course tkt4140:
        matplotlib
        numpy
        scipy
        sympy

    '''
    installed = " ... installed!"
    print ""
    print "Check for necessary modules needed for the course tkt4140"
    print ""
    
    print "check for matplotlib ",
    try:
        import matplotlib.pyplot 
        print installed
    except:
        print " IMPORT ERROR; no version of matplotlib.pyplot found"
       
        
    print "check for numpy      ",
    try:
        import numpy
        print installed
        
    except:
        print " IMPORT ERROR; no version of numpy found"
    
    print "check for scipy      ",
    try:
        import scipy
        print installed
    except:
        print " IMPORT ERROR; no version of scipy found"
        

    print "check for sympy      ",
    try:
        import sympy
        print installed
    except:
        print " IMPORT ERROR; no version of sympy found"
    

if __name__ == '__main__':
    systemCheck()
