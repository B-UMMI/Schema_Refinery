import time
import shutil
import psutil
import threading
from typing import Callable, Any

try:
    from utils import print_functions as pf
except ImportError:
    from SchemaRefinery.utils import print_functions as pf

def time_and_resource_function(monitor_memory=True, monitor_cpu=True, monitor_io=True, monitor_network=True, monitor_disk=True, monitor_threads=True, interval=0.1,) -> Callable[..., Any]:
    """
    Decorator to measure the execution time and optionally monitor resource usage.

    Parameters
    ----------
    monitor_memory : bool, optional
        Monitor memory usage, by default True
    monitor_cpu : bool, optional
        Monitor CPU usage, by default True
    monitor_io : bool, optional
        Monitor I/O operations, by default True
    monitor_network : bool, optional
        Monitor network usage, by default True
    monitor_disk : bool, optional
        Monitor disk usage, by default True
    monitor_threads : bool, optional
        Monitor number of threads, by default True
    interval : float, optional
        Monitoring interval in seconds, by default 0.1

    Returns
    -------
    Callable[..., Any]
        Wrapped function with resource monitoring.
    
    """

    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        """
        Decorator to measure the execution time and optionally monitor resource usage.

        Parameters
        ----------
        func : Callable[..., Any]
            Function to be wrapped.
        
        Returns
        -------
        Callable[..., Any]
            Wrapped function with resource monitoring.
        """
        def wrapper(*args, **kwargs) -> Any:
            process = psutil.Process()
            start_time = time.time()

            # Initialize monitoring variables
            max_memory_usage = 0
            max_cpu_usage = 0
            max_cpu_cores = 0
            max_threads = 0

            total_read_ops = 0
            total_write_ops = 0
            total_bytes_sent = 0
            total_bytes_recv = 0
            total_read_bytes = 0
            total_write_bytes = 0

            # Initial counters
            initial_io_counters = process.io_counters() if monitor_io else None
            initial_net_io_counters = psutil.net_io_counters() if monitor_network else None
            initial_disk_io_counters = psutil.disk_io_counters() if monitor_disk else None

            stop_event = threading.Event()

            def monitor_resources() -> None:
                nonlocal max_memory_usage, max_cpu_usage, max_cpu_cores, max_threads
                nonlocal total_read_ops, total_write_ops, total_bytes_sent, total_bytes_recv
                nonlocal total_read_bytes, total_write_bytes

                while not stop_event.is_set():
                    try:
                        # Monitor main process
                        total_memory_usage = process.memory_info().rss
                        total_cpu_usage = process.cpu_percent(interval=interval)
                        total_threads = process.num_threads()

                        # Monitor child processes
                        for child in process.children(recursive=True):
                            try:
                                if monitor_memory:
                                    total_memory_usage += child.memory_info().rss
                                if monitor_cpu:
                                    total_cpu_usage += child.cpu_percent(interval=0.1)
                                if monitor_threads:
                                    total_threads += child.num_threads()
                                if monitor_io:
                                    child_io_counters = child.io_counters()
                                    total_read_ops += child_io_counters.read_count
                                    total_write_ops += child_io_counters.write_count
                                if monitor_network:
                                    child_net_io_counters = child.net_io_counters()
                                    total_bytes_sent += child_net_io_counters.bytes_sent
                                    total_bytes_recv += child_net_io_counters.bytes_recv
                                if monitor_disk:
                                    try:
                                        child_disk_io_counters = child.io_counters()
                                        total_read_bytes += child_disk_io_counters.read_bytes
                                        total_write_bytes += child_disk_io_counters.write_bytes
                                    except AttributeError:
                                        pass
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                continue

                        # Update maximums
                        if monitor_memory:
                            max_memory_usage = max(max_memory_usage, total_memory_usage)
                        if monitor_cpu:
                            max_cpu_usage = max(max_cpu_usage, total_cpu_usage)
                            max_cpu_cores = max(max_cpu_cores, int(total_cpu_usage / 100 * psutil.cpu_count()))
                        if monitor_threads:
                            max_threads = max(max_threads, total_threads)

                        time.sleep(interval)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        break

            # Start resource monitoring in a separate thread
            resource_thread = threading.Thread(target=monitor_resources, daemon=True)
            resource_thread.start()

            try:
                result = func(*args, **kwargs)
            finally:
                stop_event.set()
                resource_thread.join()

            # End time
            end_time = time.time()
            execution_time = end_time - start_time

            # Print results
            terminal_width = shutil.get_terminal_size().columns
            separator = "=" * terminal_width
            pf.print_message(separator, None)
            pf.print_message("Running Stats".center(terminal_width), None)
            pf.print_message(f"Execution time: {execution_time:.2f} seconds", "info")

            if monitor_memory:
                pf.print_message(f"Maximum memory usage: {max_memory_usage / (1024 * 1024):.2f} MB", "info")
            if monitor_cpu:
                pf.print_message(f"Maximum CPU usage: {max_cpu_usage:.2f}%", "info")
                pf.print_message(f"Maximum CPU cores used: {max_cpu_cores}", "info")
                pf.print_message(f"Number of physical CPU cores: {psutil.cpu_count(logical=False)}", "info")
            if monitor_io:
                final_io_counters = process.io_counters()
                read_ops = final_io_counters.read_count - initial_io_counters.read_count + total_read_ops
                write_ops = final_io_counters.write_count - initial_io_counters.write_count + total_write_ops
                pf.print_message(f"Read operations: {read_ops}", "info")
                pf.print_message(f"Write operations: {write_ops}", "info")
            if monitor_network:
                final_net_io_counters = psutil.net_io_counters()
                bytes_sent = final_net_io_counters.bytes_sent - initial_net_io_counters.bytes_sent + total_bytes_sent
                bytes_recv = final_net_io_counters.bytes_recv - initial_net_io_counters.bytes_recv + total_bytes_recv
                pf.print_message(f"Bytes sent: {bytes_sent}", "info")
                pf.print_message(f"Bytes received: {bytes_recv}", "info")
            if monitor_disk:
                final_disk_io_counters = psutil.disk_io_counters()
                read_bytes = final_disk_io_counters.read_bytes - initial_disk_io_counters.read_bytes + total_read_bytes
                write_bytes = final_disk_io_counters.write_bytes - initial_disk_io_counters.write_bytes + total_write_bytes
                pf.print_message(f"Disk read bytes: {read_bytes}", "info")
                pf.print_message(f"Disk write bytes: {write_bytes}", "info")
            if monitor_threads:
                pf.print_message(f"Maximum number of threads: {max_threads}", "info")

            pf.print_message(separator, None)
            return result

        return wrapper

    return decorator
