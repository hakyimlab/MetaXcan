
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


class BadFilename(ReportableException):
    """Reports invalid filename"""

    def __init__(self, filename):
        super(BadFilename, self).__init__("Invalid filename: %s" % (filename))

class BadDirectory(ReportableException):
    """Reports non-existnant directory"""

    def __init__(self, dir):
        super(BadDirectory, self).__init__("Invalid directory: %s" % (dir))
