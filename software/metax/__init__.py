__version__ = "0.7.1"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
