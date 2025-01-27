import time
import shutil
import psutil
import threading
import gc
from functools import wraps
from typing import Callable, Any

try:
    from utils import print_functions as pf
except ImportError:
    from SchemaRefinery.utils import print_functions as pf

def time_and_resource_function(monitor_memory=True, monitor_cpu=True, monitor_io=True, monitor_network=True, monitor_disk=True, monitor_threads=True, monitor_gc=True, monitor_context_switches=True, monitor_open_files=True, monitor_page_faults=True, interval=0.1):
    """
    Decorator to measure the execution time and optionally the maximum memory, CPU usage, I/O operations, network usage, disk usage, number of threads, GC statistics, context switches, open files, and page faults of a function.

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
    monitor_gc : bool
        Whether to monitor garbage collection statistics. Default is True.
    monitor_context_switches : bool
        Whether to monitor context switches. Default is True.
    monitor_open_files : bool
        Whether to monitor open file descriptors. Default is True.
    monitor_page_faults : bool
        Whether to monitor page faults. Default is True.
    interval : float
        The interval in seconds between each resource measurement. Default is 0.1.

    Returns
    -------
    function
        The wrapped function with added time and optional resource measurement.
    """
    def decorator(func) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            process = psutil.Process()
            start_time = time.time()

            # Initialize monitoring variables
            max_memory_usage = 0
            max_cpu_usage = 0
            max_cpu_cores = 0
            max_threads = 0
            max_open_files = 0
            max_page_faults = 0
            total_gc_collections = [0, 0, 0]  # [generation 0, generation 1, generation 2]
            total_voluntary_context_switches = 0
            total_involuntary_context_switches = 0

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
                nonlocal max_memory_usage, max_cpu_usage, max_cpu_cores, max_threads, max_open_files, max_page_faults
                nonlocal total_gc_collections, total_voluntary_context_switches, total_involuntary_context_switches
                nonlocal total_read_ops, total_write_ops, total_bytes_sent, total_bytes_recv
                nonlocal total_read_bytes, total_write_bytes

                while not stop_event.is_set():
                    try:
                        # Monitor main process
                        total_memory_usage = process.memory_info().rss
                        total_cpu_usage = process.cpu_percent(interval=0.025) # For shorter runs use smaller interval e.g 0.025, because the children end quickly
                        total_threads = process.num_threads()
                        open_files = len(process.open_files())
                        memory_info = process.memory_full_info()
                        page_faults = getattr(memory_info, 'pfaults', 0)
                        context_switches = process.num_ctx_switches()

                        # Monitor GC statistics
                        if monitor_gc:
                            gc.collect()
                            gc_stats = gc.get_stats()
                            for i in range(3):
                                total_gc_collections[i] += gc_stats[i]['collections']

                        # Monitor child processes
                        for child in process.children(recursive=True):
                            try:
                                if monitor_memory:
                                    total_memory_usage += child.memory_info().rss
                                if monitor_cpu:
                                    total_cpu_usage += child.cpu_percent(interval=0.025) # For shorter runs use smaller interval e.g 0.025, because the children end quickly
                                if monitor_threads:
                                    total_threads += child.num_threads()
                                if monitor_io:
                                    child_io_counters = child.io_counters()
                                    total_read_ops += child_io_counters.read_count
                                    total_write_ops += child_io_counters.write_count
                                    total_read_bytes += child_io_counters.read_bytes
                                    total_write_bytes += child_io_counters.write_bytes
                                if monitor_open_files:
                                    open_files += len(child.open_files())
                                if monitor_page_faults:
                                    child_memory_info = child.memory_full_info()
                                    page_faults += getattr(child_memory_info, 'pfaults', 0)
                                if monitor_context_switches:
                                    child_context_switches = child.num_ctx_switches()
                                    total_voluntary_context_switches += child_context_switches.voluntary
                                    total_involuntary_context_switches += child_context_switches.involuntary
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                continue

                        # Update maximums
                        if monitor_memory:
                            max_memory_usage = max(max_memory_usage, total_memory_usage)
                        if monitor_cpu:
                            max_cpu_usage = max(max_cpu_usage, total_cpu_usage)
                            max_cpu_cores = max(max_cpu_cores, round(total_cpu_usage / 100))
                        if monitor_threads:
                            max_threads = max(max_threads, total_threads)
                        if monitor_open_files:
                            max_open_files = max(max_open_files, open_files)
                        if monitor_page_faults:
                            max_page_faults = max(max_page_faults, page_faults)
                        if monitor_context_switches:
                            total_voluntary_context_switches += context_switches.voluntary
                            total_involuntary_context_switches += context_switches.involuntary

                        time.sleep(interval)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue

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

            # Final counters
            final_io_counters = process.io_counters() if monitor_io else None
            final_net_io_counters = psutil.net_io_counters() if monitor_network else None
            final_disk_io_counters = psutil.disk_io_counters() if monitor_disk else None

            # Print results
            terminal_width = shutil.get_terminal_size().columns
            border_line = "=" * terminal_width
            pf.print_message(border_line, None)
            pf.print_message("Running Stats".center(terminal_width), "debug_additional_info_in_logger_only")
            pf.print_message(f"Execution time: {execution_time:.2f} seconds", "debug")

            if monitor_memory:
                pf.print_message(f"Maximum memory usage: {max_memory_usage / (1024 * 1024):.2f} MB", "debug")
            if monitor_cpu:
                pf.print_message(f"Maximum CPU usage: {max_cpu_usage:.2f}%", "debug")
                pf.print_message(f"Maximum CPU cores used: {max_cpu_cores:.2f}", "debug")
            if monitor_io and initial_io_counters and final_io_counters:
                read_ops = final_io_counters.read_count - initial_io_counters.read_count + total_read_ops
                write_ops = final_io_counters.write_count - initial_io_counters.write_count + total_write_ops
                pf.print_message(f"Read operations: {read_ops}", "debug")
                pf.print_message(f"Write operations: {write_ops}", "debug")
            if monitor_network and initial_net_io_counters and final_net_io_counters:
                # Calculate the differences for network I/O
                if initial_net_io_counters and final_net_io_counters:
                    bytes_sent = final_net_io_counters.bytes_sent - initial_net_io_counters.bytes_sent
                    bytes_recv = final_net_io_counters.bytes_recv - initial_net_io_counters.bytes_recv
                else:
                    bytes_sent = total_bytes_sent
                    bytes_recv = total_bytes_recv
                pf.print_message(f"Bytes sent: {bytes_sent}", "debug")
                pf.print_message(f"Bytes received: {bytes_recv}", "debug")
            if monitor_disk and initial_disk_io_counters and final_disk_io_counters:
                read_bytes = final_disk_io_counters.read_bytes - initial_disk_io_counters.read_bytes + total_read_bytes
                write_bytes = final_disk_io_counters.write_bytes - initial_disk_io_counters.write_bytes + total_write_bytes
                pf.print_message(f"Disk read bytes: {read_bytes}", "debug")
                pf.print_message(f"Disk write bytes: {write_bytes}", "debug")
            if monitor_threads:
                pf.print_message(f"Maximum number of threads: {max_threads}", "debug")
            if monitor_gc:
                pf.print_message(f"GC collections: {total_gc_collections}", "debug")
            if monitor_context_switches:
                pf.print_message(f"Voluntary context switches: {total_voluntary_context_switches}", "debug")
                pf.print_message(f"Involuntary context switches: {total_involuntary_context_switches}", "debug")
            if monitor_open_files:
                pf.print_message(f"Maximum open files: {max_open_files}", "debug")
            if monitor_page_faults:
                pf.print_message(f"Maximum page faults: {max_page_faults}", "debug")

            pf.print_message(border_line, "debug_additional_info_in_logger_only")
            return result

        return wrapper

    return decorator

def time_function():
    """
    Decorator to measure the execution time of a function.

    Parameters
    ----------
    None

    Returns
    -------
    function
        The wrapped function with added time measurement.
    """
    def decorator(func) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> Any:
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time
            pf.print_message(f"Execution time: {execution_time:.2f} seconds", "info")
            return result

        return wrapper

    return decorator