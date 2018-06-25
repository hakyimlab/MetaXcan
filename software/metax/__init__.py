__version__ = "0.6.4"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
