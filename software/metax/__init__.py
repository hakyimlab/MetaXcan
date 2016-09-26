__version__ = "0.3.2"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
