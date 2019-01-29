__version__ = "0.6.8"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
