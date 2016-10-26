__version__ = "0.3.3"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
