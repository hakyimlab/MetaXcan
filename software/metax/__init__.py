__version__ = "0.7.4"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
