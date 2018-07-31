__version__ = "0.6.5"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
