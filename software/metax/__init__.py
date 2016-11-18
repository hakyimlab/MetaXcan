__version__ = "0.3.4"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
