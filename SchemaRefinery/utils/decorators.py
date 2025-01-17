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

def time_and_resource_function(monitor_memory=True, monitor_cpu=True, monitor_io=True):
    """
    Decorator to measure the execution time and optionally the maximum memory, CPU usage, and I/O operations of a function.

    Parameters
    ----------
    monitor_memory : bool
        Whether to monitor memory usage. Default is True.
    monitor_cpu : bool
        Whether to monitor CPU usage. Default is True.
    monitor_io : bool
        Whether to monitor I/O operations. Default is True.

    Returns
    -------
    function
        The wrapped function with added time and optional resource measurement.
    """
    def decorator(func) -> Callable[..., Any]:
        """
        The decorator function that measures the execution time and resources.

        Parameters
        ----------
        func : function
            The function to be decorated.
        
        Returns
        -------
        function
            The wrapped function that measures the execution time and
        """
        def wrapper(*args, **kwargs) -> Any:
            """
            The wrapper function that measures the execution time and resources.

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
            max_cpu_usage = 0
            initial_io_counters = process.io_counters() if monitor_io else None
            stop_event = threading.Event()

            def monitor_resources() -> None:
                """
                Monitor the memory and CPU usage of the process.

                This function is run in a separate thread to monitor the memory and CPU usage of the process.
                
                Returns
                -------
                None
                """
                nonlocal max_memory_usage, max_cpu_usage
                while not stop_event.is_set():
                    if monitor_memory:
                        mem_info = process.memory_info()
                        max_memory_usage = max(max_memory_usage, mem_info.rss)
                    if monitor_cpu:
                        cpu_usage = process.cpu_percent(interval=0.1)
                        max_cpu_usage = max(max_cpu_usage, cpu_usage)
                    time.sleep(0.1)

            resource_thread = threading.Thread(target=monitor_resources)
            resource_thread.start()

            try:
                result = func(*args, **kwargs)
            finally:
                stop_event.set()
                resource_thread.join()

            end_time = time.time()
            execution_time = end_time - start_time
            pf.print_message(f"Execution time: {execution_time:.2f} seconds", "info")
            if monitor_memory:
                pf.print_message(f"Maximum memory usage: {max_memory_usage / (1024 * 1024):.2f} MB", "info")
            if monitor_cpu:
                pf.print_message(f"Maximum CPU usage: {max_cpu_usage:.2f}%", "info")
            if monitor_io:
                final_io_counters = process.io_counters()
                read_ops = final_io_counters.read_count - initial_io_counters.read_count
                write_ops = final_io_counters.write_count - initial_io_counters.write_count
                pf.print_message(f"Read operations: {read_ops}", "info")
                pf.print_message(f"Write operations: {write_ops}", "info")
            return result

        return wrapper
    return decorator