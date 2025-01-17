import time
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
