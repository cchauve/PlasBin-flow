import sys
import os
import logging
import subprocess

def check_file(in_file):
    if not os.path.isfile(in_file):
        msg = f'{in_file} is missing'
        logging.error(msg)
        print(f'ERROR\t{msg}', file=sys.stderr)
        sys.exit(1)

def log_file(in_file):
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
        logging.exception(f'Running {cmd_str}')
        print(f'ERROR\tRunning {cmd_str}', file=sys.stderr)
        sys.exit(1)
    else:
        logging.info(f'STDOUT:\n{process.stdout}')
        logging.warning(f'STDERR:\n{process.stderr}')      
