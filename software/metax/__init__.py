__version__ = "0.5.9"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
