__version__ = "0.3.1"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
