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

def run_cmd(cmd):
    """ Run external command """
    cmd_str = ' '.join(cmd)
    logging.info(f'COMMAND {cmd_str}')
    try:
        process = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        msg = f'Running {cmd_str}: {e}'
        process_exception(msg)
    else:
        logging.info(f'STDOUT:\n{process.stdout}')
        logging.warning(f'STDERR:\n{process.stderr}')      
