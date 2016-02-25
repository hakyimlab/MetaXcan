#! /usr/bin/env python
__author__ = 'heroico'

import Tkinter
import ttk
import metax.MainScreen as MainScreen
import metax.Logging as Logging
class UI(object):
    def __init__(self, root):
        self.controller = MainScreen.MainScreen(root, self)

#shamelessly stolen from the interwebz: http://stackoverflow.com/questions/3352918/how-to-center-a-window-on-the-screen-in-tkinter
def center(toplevel):
    toplevel.update_idletasks()
    w = toplevel.winfo_screenwidth()
    h = toplevel.winfo_screenheight()
    size = tuple(int(_) for _ in toplevel.geometry().split('+')[0].split('x'))
    x = w/2 - size[0]/2
    y = h/2 - size[1]/2
    toplevel.geometry("%dx%d+%d+%d" % (size + (x, y)))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='UI for MetaXcan.')

    parser.add_argument("--verbosity",
                        help="Log verbosity level. 1 is everything being logged. 10 is only high level messages, above 10 will hardly log anything",
                        default = "10")

    args = parser.parse_args()
    Logging.configureLogging(int(args.verbosity))

    root = Tkinter.Tk()
    style = ttk.Style()
    style.theme_use('clam') #clam, default, alt
    #root.resizable(width=False, height=False)
    #root.geometry('800x600')
    app = UI(root)
    center(root)
    root.mainloop()
    #root.destroy()



