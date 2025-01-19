import logging

def setup_logger(log_file):
    """
    Set up a logger object with a file handler.

    Parameters
    ----------
    log_file : str
        The path to the log file.
    
    Returns
    -------
    logging.Logger
        The logger object.
    """
    logger = logging.getLogger('SchemaRefinery')
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger

def get_log_file_path(logger):
    """
    Retrieve the path to the log file from the logger object.

    Parameters
    ----------
    logger : logging.Logger
        The logger object from which to retrieve the log file path.

    Returns
    -------
    str
        The path to the log file, or None if no file handler is found.
    """
    if logger is None:
        return None
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            return handler.baseFilename
    return None