import os

def create_directory(dir:str):
    if not os.path.isdir(dir):
        os.mkdir(dir)

def check_and_delete_file(file:str):
    if os.path.isfile(file):
        os.remove(file)
