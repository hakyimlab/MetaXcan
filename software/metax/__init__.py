__version__ = "0.7.2"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
