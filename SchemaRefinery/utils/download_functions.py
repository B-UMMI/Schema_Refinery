import urllib.request
import time
import http.client
from typing import Union, Tuple

def download_file(url: str, file_name: str, retry: int) -> Union[str, Tuple[str, Any]]:
    """
    Downloads a file from the given URL and saves it with the specified file name. Retries the download
    up to the specified number of times if it fails.

    Parameters
    ----------
    url : str
        The URL to download the file from.
    file_name : str
        The name of the file to be saved.
    retry : int
        Maximum number of retries if the download fails.

    Returns
    -------
    Union[str, Tuple[str, Any]]
        A string indicating that the download failed or a tuple with the response information for a successful download.
    """

    tries: int = 0
    response: Union[str, Tuple[str, http.client.HTTPMessage]] = ('', None)
    while tries < retry:
        try:
            response = urllib.request.urlretrieve(url, file_name)
            break
        except Exception:
            response = f'Failed: {file_name}'
            tries += 1
            print(f'Retrying {file_name.split("/")[-1]} ...{tries}')
            time.sleep(1)

    return response
