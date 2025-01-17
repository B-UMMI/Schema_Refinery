import time
import shutil
import psutil
import threading
from typing import Callable, Any

try:
    from utils import print_functions as pf
except ImportError:
    from SchemaRefinery.utils import print_functions as pf

def time_and_resource_function(monitor_memory=True, monitor_cpu=True, monitor_io=True, monitor_network=True, monitor_disk=True, monitor_threads=True):
    """
    Decorator to measure the execution time and optionally the maximum memory, CPU usage, I/O operations, network usage, disk usage, and number of threads of a function.

    Parameters
    ----------
    monitor_memory : bool
        Whether to monitor memory usage. Default is True.
    monitor_cpu : bool
        Whether to monitor CPU usage. Default is True.
    monitor_io : bool
        Whether to monitor I/O operations. Default is True.
    monitor_network : bool
        Whether to monitor network usage. Default is True.
    monitor_disk : bool
        Whether to monitor disk usage. Default is True.
    monitor_threads : bool
        Whether to monitor the number of threads. Default is True.

    Returns
    -------
    function
        The wrapped function with added time and optional resource measurement.
    """
    def decorator(func) -> Callable[..., Any]:
        def wrapper(*args, **kwargs) -> Any:
            process = psutil.Process()
            start_time = time.time()
            max_memory_usage = 0
            max_cpu_usage = 0
            initial_io_counters = process.io_counters() if monitor_io else None
            initial_net_io_counters = psutil.net_io_counters() if monitor_network else None
            initial_disk_io_counters = psutil.disk_io_counters() if monitor_disk else None
            max_threads = 0
            stop_event = threading.Event()

            def monitor_resources() -> None:
                nonlocal max_memory_usage, max_cpu_usage, max_threads
                while not stop_event.is_set():
                    if monitor_memory:
                        mem_info = process.memory_info()
                        max_memory_usage = max(max_memory_usage, mem_info.rss)
                    if monitor_cpu:
                        cpu_usage = process.cpu_percent(interval=0.1)
                        max_cpu_usage = max(max_cpu_usage, cpu_usage)
                    if monitor_threads:
                        max_threads = max(max_threads, process.num_threads())
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

            terminal_width = shutil.get_terminal_size().columns
            separator = "=" * terminal_width
            pf.print_message(separator, None)
            pf.print_message("Running Stats".center(terminal_width), None)
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
            if monitor_network:
                final_net_io_counters = psutil.net_io_counters()
                bytes_sent = final_net_io_counters.bytes_sent - initial_net_io_counters.bytes_sent
                bytes_recv = final_net_io_counters.bytes_recv - initial_net_io_counters.bytes_recv
                pf.print_message(f"Bytes sent: {bytes_sent}", "info")
                pf.print_message(f"Bytes received: {bytes_recv}", "info")
            if monitor_disk:
                final_disk_io_counters = psutil.disk_io_counters()
                read_bytes = final_disk_io_counters.read_bytes - initial_disk_io_counters.read_bytes
                write_bytes = final_disk_io_counters.write_bytes - initial_disk_io_counters.write_bytes
                pf.print_message(f"Disk read bytes: {read_bytes}", "info")
                pf.print_message(f"Disk write bytes: {write_bytes}", "info")
            if monitor_threads:
                pf.print_message(f"Maximum number of threads: {max_threads}", "info")
            pf.print_message(separator, None)
            return result

        return wrapper
    return decorator