import shutil
from datetime import datetime

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

def print_message(message: str, message_type: str = "info", logger: bool = False, end = '\n', flush: bool = False) -> None:
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

    print(formatted_message, end = end, flush = flush)
