__version__ = "0.7.5"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
