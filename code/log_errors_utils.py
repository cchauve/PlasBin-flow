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
        
def _check_file(in_file, log=False, msg='FILE'):
    try:
        if not os.path.isfile(in_file):
            raise CustomException('File is missing')
        elif os.path.getsize(in_file) == 0:
            process_warning(f'{msg}\t{in_file}: is empty')
    except Exception as e:
        process_exception(f'{msg}\t{in_file}: {e}')
    else:
        if log:
            logging.info(f'{msg}\t{in_file}')
            
def check_file(in_file):
    _check_file(in_file, log=False)

def log_file(in_file):
    _check_file(in_file, log=True)
        
def check_number_range(x, allowed_range=(None,None), msg='CHECK_RANGE'):
    """
    Check that number x is in an allowed range
    Args:
        - x: int or float
        - allowed_range = List((None or number),(None or number))
          if a boundary of the range is None, it is not checked
    """
    try:
        y,z = allowed_range[0],allowed_range[1]
        if y is not None and x < y:
            raise CustomException(
                f'{x} outside of allowed range: {x} < {y}'
            )
        if z is not None and x > z:
            raise CustomException(
                f'{x} outside of allowed range: {x} > {z}'
            )
    except Exception as e:
        process_exception(f'{msg}: {e}')
        
def check_number_eq(x, y, msg='CHECK_EQ'):
    try:
        if x != y:
            raise CustomException(f'{x} != {y}')
    except Exception as e:
        process_exception(f'{msg}: {e}')

def check_number_gt(x, y, msg='CHECK_EQ'):
    try:
        if x <= y:
            raise CustomException(f'{x} <= {y}')
    except Exception as e:
        process_exception(f'{msg}: {e}')

def check_number_lt(x, y, msg='CHECK_EQ'):
    try:
        if x >= y:
            raise CustomException(f'{x} >= {y}')
    except Exception as e:
        process_exception(f'{msg}: {e}')    
        
def check_lists_equality(list1, list2, msg='CHECK_LISTS_EQ'):
    try:
        if set(list1) != set(list2):
            raise CustomException(f'Sets equality error')
    except Exception as e:
        process_exception(f'{msg}: {e}')

def check_lists_inclusion(list1, list2, msg='CHECK_LISTS_INCLUSION'):
    try:
        if not set(list1).issubset(list2):
            raise CustomException(f'Sets inclusion error')
    except Exception as e:
        process_exception(f'{msg}: {e}')    

def check_num_fields(in_list, min_num_fields, msg='CHECK_NUM_FIELDS'):
    try:
        num_fields = len(in_list)
        if num_fields < min_num_fields:
            raise CustomException(
                f'At least {min_num_fields} expected fields, {num_fields} read'
            )
    except Exception as e:
        process_exception(f'{msg}: {e}')
        
def check_num_fields_eq(in_list, exact_num_fields, msg='CHECK_NUM_FIELDS'):
    try:
        num_fields = len(in_list)
        if num_fields != exact_num_fields:
            raise CustomException(
                f'Eaxctly {min_num_fields} expected fields, {num_fields} read'
            )
    except Exception as e:
        process_exception(f'{msg}: {e}')
    
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

""" Files and directories function """
        
def clean_files(files2clean):
    for in_file in files2clean:
        if os.path.isfile(in_file):
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

