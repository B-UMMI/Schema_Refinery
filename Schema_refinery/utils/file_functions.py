import os

def check_and_make_directory(dir:str):
    if not os.path.isdir(dir):
        os.mkdir(dir)