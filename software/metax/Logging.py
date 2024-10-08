__author__ = 'heroico'

import logging
import sys

def configureLogging(level=5, target=sys.stderr):
    logger = logging.getLogger()
    logger.setLevel(level)

    kludge=False
    if logger.hasHandlers():
        kludge=True
        logger.handlers.clear()

    # create console handler and set level to info
    RED = "\033[91m"
    RESET = "\033[0m"

    # A custom logging formatter for warning
    class CustomFormatter(logging.Formatter):
        def format(self, record):
            if record.levelno == logging.WARNING:
                record.msg = f"{RED}{record.msg}{RESET}"
            return super().format(record)

    # create console handler and set level to info
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = CustomFormatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if kludge:
        logging.info("Had to reset logging handlers. Why?")