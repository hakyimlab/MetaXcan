
def silentRm(filename):
    try:
        os.unlink(filename)
    except:
        pass
