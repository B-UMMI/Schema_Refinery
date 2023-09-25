import os

def create_directory(dir:str):
    if not os.path.isdir(dir):
        os.mkdir(dir)
