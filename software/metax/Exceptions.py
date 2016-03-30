
class ReportableException(Exception):
    """Simple exeception with message"""
    def __init__(self, msg):
        self.msg = msg

class InvalidInputFormat(ReportableException):
    """Error associated with input format"""
    def __init__(self, msg):
        super(InvalidInputFormat, self).__init__("Invalid input format: %s" % (msg))

class InvalidOutputFormat(ReportableException):
    """Error associated with input format"""
    def __init__(self, msg):
        super(InvalidOutputFormat, self).__init__("Invalid output format: %s" % (msg))

class MalformedInputFile(ReportableException):
    """Errors associated with a problematic input file"""
    def __init__(self, filename, msg):
        super(MalformedInputFile, self).__init__("While processing file, %s: %s" % (filename, msg) )

class BadFilename(ReportableException):
    """Reports invalid filename"""

    def __init__(self, filename):
        super(BadFilename, self).__init__("Invalid filename: %s" % (filename))

class BadDirectory(ReportableException):
    """Reports non-existnant directory"""

    def __init__(self, dir):
        super(BadDirectory, self).__init__("Invalid directory: %s" % (dir))

class InvalidDbFormat(ReportableException):
    """Report invalid table names, missing columns, etc"""

    def __init__(self, filename, msg):
        errorMsg = "DB Error encounted loading data from %s\n\t'%s'" % (filename, msg)
        super(InvalidDbFormat, self).__init__(errorMsg)