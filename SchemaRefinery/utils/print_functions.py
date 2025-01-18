import os
import time
import shutil
import psutil
import logging
import platform
import subprocess
from importlib.metadata import version, PackageNotFoundError
from datetime import datetime
from typing import Optional

try:
    from utils import (globals as gb,
                       constants as ct)
except ModuleNotFoundError:
    from SchemaRefinery.utils import (globals as gb,
                                      constants as ct)

def print_logo() -> None:
    """
    Print the SchemaRefinery logo.

    This function prints the SchemaRefinery logo in the terminal.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    big_logo = """


███████╗ ██████╗██╗  ██╗███████╗███╗   ███╗ █████╗ ██████╗ ███████╗███████╗██╗███╗   ██╗███████╗██████╗ ██╗   ██╗
██╔════╝██╔════╝██║  ██║██╔════╝████╗ ████║██╔══██╗██╔══██╗██╔════╝██╔════╝██║████╗  ██║██╔════╝██╔══██╗╚██╗ ██╔╝
███████╗██║     ███████║█████╗  ██╔████╔██║███████║██████╔╝█████╗  █████╗  ██║██╔██╗ ██║█████╗  ██████╔╝ ╚████╔╝ 
╚════██║██║     ██╔══██║██╔══╝  ██║╚██╔╝██║██╔══██║██╔══██╗██╔══╝  ██╔══╝  ██║██║╚██╗██║██╔══╝  ██╔══██╗  ╚██╔╝  
███████║╚██████╗██║  ██║███████╗██║ ╚═╝ ██║██║  ██║██║  ██║███████╗██║     ██║██║ ╚████║███████╗██║  ██║   ██║   
╚══════╝ ╚═════╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝   ╚═╝   
                                                                                                                 
                                                                                                                 

    """

    small_logo = """

███████╗██████╗ 
██╔════╝██╔══██╗
███████╗██████╔╝
╚════██║██╔══██╗
███████║██║  ██║
╚══════╝╚═╝  ╚═╝
                
    """

    size = shutil.get_terminal_size()
    if size.columns < 100:
        # Center the small logo
        for line in small_logo.splitlines():
            print(line.center(size.columns))
    else:
        # Center the big logo
        for line in big_logo.splitlines():
            print(line.center(size.columns))

def print_module_currently_running(module: str) -> None:
    """
    Print the module that is currently running.

    This function prints the name of the module that is currently running.

    Parameters
    ----------
    module : str
        The name of the module that is currently running.

    Returns
    -------
    None
    """
    size = shutil.get_terminal_size()
    box_width = size.columns
    module_line = f"{module}".center(box_width)
    border_line = "=" * box_width

    to_print = f"{border_line}\n{module_line}\n{border_line}\n"

    print(to_print)

def print_message(message: str, message_type: Optional[str] = "info", end = '\n', flush: bool = False) -> None:
    """
    Print a formatted message with the current time and message type.

    Parameters
    ----------
    message : str
        The message to print.
    message_type : str
        The type of message (e.g., "info", "warning", "error").
    logger : bool
        Whether to log the message.

    Returns
    -------
    None
    """
    logger: Optional[logging.Logger] = gb.LOGGER # This can change to a logger object if user uses --logger

    current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if message_type == "info":
        formatted_message = f"[INFO] {current_time} - {message}"
    elif message_type == "warning":
        formatted_message = f"[WARNING] {current_time} - {message}"
    elif message_type == "error":
        formatted_message = f"[ERROR] {current_time} - {message}"
    elif message_type == "debug":
        formatted_message = f"[DEBUG] {current_time} - {message}"
    elif message_type == None:
        formatted_message = f"{message}"
    else:
        formatted_message = f"{current_time} - {message}"
    # Log the message if logger is not None
    if logger:
        if message_type == "info":
            logger.info(message)
        elif message_type == "warning":
            logger.warning(message)
        elif message_type == "error":
            logger.error(message)
        elif message_type == "debug":
            logger.debug(message)
        else:
            logger.info(message)
    
    print(formatted_message, end = end, flush = flush)

def print_system_info() -> None:
    """
    Print the system information.

    This function prints the system information.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    def get_processor_name():
        try:
            result = subprocess.run(["lscpu"], capture_output=True, text=True, check=True)
            for line in result.stdout.splitlines():
                if "Model name" in line:
                    return line.split(":")[1].strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return platform.processor()

    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print_message(separator, None)
    print_message("System Information".center(terminal_width), None)
    print_message(f"Operating System: {platform.system()} {platform.release()}", 'info')
    print_message(f"OS Version: {platform.version()}", 'info')
    print_message(f"Machine: {platform.machine()}", 'info')
    print_message(f"Processor: {get_processor_name()}", 'info')
    print_message(f"CPU count: {psutil.cpu_count(logical=True)} (logical), {psutil.cpu_count(logical=False)} (physical)", 'info')
    print_message(f"Total memory: {psutil.virtual_memory().total / (1024 * 1024):.2f} MB", 'info')
    print_message(f"Available memory: {psutil.virtual_memory().available / (1024 * 1024):.2f} MB", 'info')
    print_message(f"Disk usage: {psutil.disk_usage('/').percent}%", 'info')
    # Calculate and print system uptime in days, hours, and seconds
    current_time = time.time()
    boot_time = psutil.boot_time()
    uptime_seconds = current_time - boot_time

    uptime_days = int(uptime_seconds // (24 * 3600))
    uptime_hours = int((uptime_seconds % (24 * 3600)) // 3600)
    uptime_minutes = int((uptime_seconds % 3600) // 60)
    uptime_remaining_seconds = int(uptime_seconds % 60)
    print_message(f"System Uptime: {uptime_days} days, {uptime_hours} hours, {uptime_minutes} minutes, {uptime_remaining_seconds} seconds", 'info')
    print_message(separator, None)

def print_schema_refinery_info():
    """
    Print the SchemaRefinery information.

    This function prints the SchemaRefinery information.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    schema_refinery_path = os.path.dirname(os.path.abspath(__file__))

    print_message(separator, None)
    print_message("SchemaRefinery Information".center(terminal_width), None)
    print_message(f"SchemaRefinery Version: {ct.VERSION}", 'info')
    print_message(f"SchemaRefinery Path: {schema_refinery_path}", 'info')
    print_message(separator, None)

def print_dependencies_info(dependencies):
    """
    Print the dependencies information.

    This function prints the dependencies information.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    terminal_width = shutil.get_terminal_size().columns
    separator = "=" * terminal_width

    print_message(separator, None)
    print_message("Dependencies Information".center(terminal_width), None)
    print_message(f"Python Version: {platform.python_version()}", 'info')
    for dep in dependencies:
        try:
            dep_version = version(dep)
            print_message(f"{dep} version: {dep_version}", 'info')
        except PackageNotFoundError:
            print_message(f"{dep} is not installed", 'warning')

    # Check if BLAST is installed
    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True)
        if result.returncode == 0:
            print_message(f"BLAST is installed: {result.stdout.splitlines()[0]}", 'info')
        else:
            print_message("BLAST is not installed", 'warning')
    except FileNotFoundError:
        print_message("BLAST is not installed", 'warning')

    # Check if Datasets is installed
    try:
        result = subprocess.run(["datasets", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            print_message(f"Datasets is installed: {result.stdout.strip()}", 'info')
        else:
            print_message("Datasets is not installed", 'warning')
    except FileNotFoundError:
        print_message("Datasets is not installed", 'warning')

    print_message(separator, None)