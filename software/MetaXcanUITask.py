__author__ = 'heroico'

# most of this stolen from http://stackoverflow.com/questions/3842155/is-there-a-way-to-make-the-tkinter-text-widget-read-only
import time
import logging
from multiprocessing import Process
from multiprocessing.queues import Queue
from threading import Thread
from Tkinter import *
import Logging

class Monitor(Thread):
    def __init__(self, callback, queue):
        super(Monitor, self).__init__(group=None, target=None, name=None, args=(), kwargs={})
        self.callback = callback
        self.queue = queue
        self.give_up = False
        self.daemon = True

    def run(self):
        while not self.give_up:
            self.callback(self.queue)

# This is a Queue that behaves like stdout
class StdQueue(Queue):
    def __init__(self,*args,**kwargs):
        Queue.__init__(self,*args,**kwargs)

    def write(self,msg):
        self.put(msg)

    def flush(self):
        sys.__stdout__.flush()

def _dispatch(q, work):
    sys.stdout = q
    Logging.configureLogging(10, sys.stdout)
    work.run()

def runWork(work, callback):
    q = StdQueue()

    monitor = Monitor(callback, q)
    monitor.start()

    process = Process(target=_dispatch, args=(q,work,))
    process.start()

    return process, monitor

class Dummy(object):
    def __init__(self,i):
        self.i = i

    def run(self):
        logging.info("running")
        time.sleep(3)
        # for i in xrange(self.i,self.i+1500000):
        #     logging.info("%i", i)
        logging.info("ran")
        exit(0)
