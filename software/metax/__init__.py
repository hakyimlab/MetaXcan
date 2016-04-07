__version__ = "0.2.3"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
