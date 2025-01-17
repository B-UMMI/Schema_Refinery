import time
import psutil
import threading
from typing import Callable, Any

try:
    from utils import print_functions as pf
except ImportError:
    from SchemaRefinery.utils import print_functions as pf

def time_function(func: Callable[..., Any]) -> Callable[..., Any]:
    """
    A decorator that measures and prints the execution time of a function.

    Parameters
    ----------
    func : Callable[..., Any]
        The function to be timed.

    Returns
    -------
    Callable[..., Any]
        The wrapped function that measures and prints its execution time.
    """
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        """
        The wrapper function that measures the execution time.

        Parameters
        ----------
        *args : Any
            Positional arguments for the decorated function.
        **kwargs : Any
            Keyword arguments for the decorated function.

        Returns
        -------
        Any
            The result of the decorated function.
        """
        # Record the start time
        start_time = time.time()
        # Call the original function and get the result
        result = func(*args, **kwargs)
        # Record the end time
        end_time = time.time()
        pf.print_message(f"Execution time: {end_time - start_time:.2f} seconds", message_type="info")
        return result
    return wrapper

def time_and_memory_function(func) -> Callable[..., Any]:
    """
    A decorator that measures and prints the execution time and memory usage of a function.

    Parameters
    ----------
    func : Callable[..., Any]
        The function to be timed.

    Returns
    -------
    Callable[..., Any]
        The wrapped function that measures and prints its execution time and memory usage
    """
    def wrapper(*args, **kwargs) -> Any:
        """
        The wrapper function that measures the execution time and memory usage.

        Parameters
        ----------
        *args : Any
            Positional arguments for the decorated function.
        **kwargs : Any
            Keyword arguments for the decorated function.
        
        Returns
        -------
        Any
            The result of the decorated function.
        """
        process = psutil.Process()
        start_time = time.time()
        max_memory_usage = 0
        stop_event = threading.Event()

        def monitor_memory() -> None:
            """
            Monitor the memory usage of the process.

            Returns
            -------
            None
            """
            nonlocal max_memory_usage
            while not stop_event.is_set():
                mem_info = process.memory_info()
                max_memory_usage = max(max_memory_usage, mem_info.rss)
                time.sleep(0.1)

        memory_thread = threading.Thread(target=monitor_memory)
        memory_thread.start()

        try:
            result = func(*args, **kwargs)
        finally:
            stop_event.set()
            memory_thread.join()

        end_time = time.time()
        execution_time = end_time - start_time
        pf.print_message(f"Execution time: {execution_time:.2f} seconds", "info")
        pf.print_message(f"Maximum memory usage: {max_memory_usage / (1024 * 1024):.2f} MB", "info")
        return result

    return wrapper