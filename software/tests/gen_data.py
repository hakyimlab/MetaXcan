
import sys

class StructuredExample(object):
    def __init__(self, prefix="__structured_test", header=None, dtypes=None):
        self.filename = "%s.txt" % (prefix)
        self.header = header

        if header is not None and dtypes is not None:
            if len(header) != len(dtypes):
                print("Header and dtypes must be the same size: ", header, dtypes)
                sys.exit(1)


#def GenerateStructuredExample(filename="__file_struct")