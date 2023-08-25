import sys
import os
import logging

def _check_file(in_file):
    if not os.path.isfile(in_file):
        msg = f'{in_file} is missing'
        logging.error(msg)
        print(f'ERROR\t{msg}', file=sys.stderr)
        sys.exit(1)

def _log_file(in_file):
    """ Write logging message for creating file in_file """
    if os.path.isfile(in_file):
        logging.info(f'FILE\t{in_file}')
    else:
        logging.error(f'FILE\t{in_file} is missing')
        print(f'ERROR\t{in_file} is missing', file=sys.stderr)
        sys.exit(1)

class CustomException(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)
