import sys
import os
import logging
import subprocess

class CustomException(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)

def process_exception(msg):
    logging.exception(msg)
    print(f'ERROR\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_error(msg):
    logging.error(msg)
    print(f'ERROR\t{msg}', file=sys.stderr)
    sys.exit(1)

def check_file(in_file):
    if not os.path.isfile(in_file):
        process_error(f'{in_file} is missing')

def log_file(in_file):
    """ Write logging message for creating file in_file """
    if os.path.isfile(in_file):
        logging.info(f'FILE\t{in_file}')
    else:
        process_error(f'File\t{in_file} is missing')

def clean_files(files2clean):
    for in_file in files2clean:
        if os.path.isfile(in_file):
            os.remove(in_file)
        
def create_directory(in_dir_list):
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

def _run_cmd(cmd, output, num_attempts, exit_on_error=True):
    """ 
    Run external command, trying at most num_attempts times  
    If output is None, write output in logging file
    """
    cmd_str = ' '.join(cmd)
    logging.info(f'COMMAND {cmd_str}')
    attempt = 1
    process_returncode = -1
    while attempt <= num_attempts:
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            msg = f'Running {cmd_str}: {e} attempt #{attempt}'
            if attempt < num_attempts:
                logging.warning(f'{msg}: retrying')
            elif exit_on_error:
                process_exception(f'{msg}: aborting')
            else:
                logging.error(f'{msg}: failed but not aborting')
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
    
def run_cmd(cmd, num_attempts=5):
    """ Run external command, trying at most num_attempts=5 times  """
    return _run_cmd(cmd, None, num_attempts)
    
def run_cmd_redirect(cmd, out_file_name, num_attempts=5):
    """ Run external command with redirection, trying at most num_attempts=5 times  """
    return _run_cmd(cmd, out_file_name, num_attempts)

