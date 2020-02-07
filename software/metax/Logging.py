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
    handler = logging.StreamHandler(target)
    handler.setLevel(level)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if kludge:
        logging.info("Had to reset logging handlers. Why?")