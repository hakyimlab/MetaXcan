__version__ = "0.6.12"

def exitIf(doExit, Exception, msg):
    if doExit:
        raise Exception(msg)
