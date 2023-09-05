""" Handling logging, errors and running external commands """

import sys
import os
import logging
import subprocess

# Maximum number of attempts for an external command
DEFAULT_RUN_ATTEMPTS = 5

""" Exceptions handling """

class CustomException(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)

EXCEPTION_MISSING_FILE = CustomException('File is missing')
EXCEPTION_EMPTY_FILE = CustomException('File is empty')
EXCEPTION_BELOW_RANGE = CustomException('Lower than allowed range')
EXCEPTION_ABOVE_RANGE = CustomException('Greater than allowed range')
EXCEPTION_UNEQUAL_SETS = CustomException('Inconsistents sets')
EXCEPTION_EMPTY_SET = CustomException('Empty set')
        
def check_file(in_file):
    if not os.path.isfile(in_file):
        raise EXCEPTION_MISSING_FILE
    elif os.path.getsize(in_file) == 0:
        raise EXCEPTION_EMPTY_FILE

def check_number(x, allowed_range=(0,1)):
    """
    Check that number x is in an allowed range
    Args:
        - x: int or float
        - allowed_range = List((None or number),(None or number))
          if a boundary of the range is None, it is not checked
    """
    y,z = allowed_range[0],allowed_range[1]
    if y is not None and x < y:
        raise EXCEPTION_BELOW_RANGE
    if z is not None and x > z:
        raise EXCEPTION_ABOVE_RANGE

def check_lists(list1, list2):
    if sorted(list1) != sorted(list2):
        raise EXCEPTION_UNEQUAL_SETS

def process_exception(msg):
    logging.exception(msg)
    print(f'EXCEPTION\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_error(msg):
    logging.error(msg)
    print(f'ERROR\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_warning(msg):
    logging.warning(msg)
    print(f'WARNING\t{msg}', file=sys.stderr)
    
""" Logging functions """
        
def log_file(in_file):
    """ Write logging message for creating file in_file """
    try:
        check_file(in_file)
    except EXCEPTION_MISSING_FILE as e:
        process_exception(f'FILE\t{in_file}: {e}')
    except EXCEPTION_EMPTY_FILE as e:
        process_warning(f'FILE\t{in_file}: {e}')
    else:
        logging.info(f'FILE\t{in_file}')

""" Files and directories function """
        
def clean_files(files2clean, msg='Deleting file'):
    for in_file in files2clean:
        try:
            check_file(in_file)
        except EXCEPTION_MISSING_FILE as e:
            process_exception(f'{msg}, {in_file}: {e}')
        else:
            os.remove(in_file)

def create_directory(in_dir_list):
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

""" Subprocess functions """
            
def _run_cmd(cmd, output, num_attempts, exit_on_error):
    """ 
    Run external command, trying at most num_attempts times  
    If output is None, write output in logging file
    """
    cmd_str = ' '.join(cmd)
    logging.info(f'COMMAND\t{cmd_str}')
    attempt = 1
    process_returncode = -1
    while attempt <= num_attempts:
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            msg = f'COMMAND\t{cmd_str} attempt #{attempt} {e}'
            if attempt < num_attempts:
                process_warning(f'{msg}: retrying')
            elif exit_on_error:
                process_exception(f'{msg}: aborting')
            else:
                process_warning(f'{msg}: failed but not aborting')
        else:
            if output is None:
                logging.info(f'STDOUT:\n{process.stdout}')
            else:
                with open(output, 'w') as out_file:                    
                    out_file.write(process.stdout)
            if len(process.stderr) > 0:
                logging.warning(f'STDERR:\n{process.stderr}')
            process_returncode = process.returncode
        attempt += 1
    return process_returncode
    
def run_cmd(cmd, num_attempts=DEFAULT_RUN_ATTEMPTS, exit_on_error=True):
    """ Run external command, trying at most num_attempts=5 times  """
    return _run_cmd(cmd, None, num_attempts, exit_on_error)
    
def run_cmd_redirect(cmd, out_file_name, num_attempts=DEFAULT_RUN_ATTEMPTS, exit_on_error=True):
    """ Run external command with redirection, trying at most num_attempts=5 times  """
    return _run_cmd(cmd, out_file_name, num_attempts, exit_on_error)

