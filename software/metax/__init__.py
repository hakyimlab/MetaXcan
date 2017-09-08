__version__ = "0.5.6"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
