__version__ = "0.2.5"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
