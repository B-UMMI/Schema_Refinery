import urllib.request
import time

def download_file(url, file_name, retry):
    """Accept a URL to download a file.

    Parameters
    ----------
    url : str
        An url to download a file.
    file_name : str
        The name of the file to be downloaded.
    retry : int
        Maximum number of retries if download fails.

    Returns
    -------
    response : str
        A string indicating that the download failed or
        an object with the response information for a
        successful download.
    """

    tries = 0
    while tries < retry:
        try:
            response = urllib.request.urlretrieve(url, file_name)
            break
        except Exception:
            response = 'Failed: {0}'.format(file_name)
            tries += 1
            print('Retrying {0} ...{1}'.format(file_name.split('/')[-1], tries))
            time.sleep(1)

    return response
