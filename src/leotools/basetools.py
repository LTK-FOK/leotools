""""This module contains tools supplmentary to geoinformatic applications."""

from pathlib import Path
import time
import fnmatch
import shutil
import zipfile
import tarfile

### Basic stuff

def timestamp():
    """Returns the timestamp as a string."""
    time_format = '%Y-%m-%d %H:%M:%S'
    ts = time.strftime(time_format, time.localtime(time.time()))
    return ts

def format_time(a):
    """Formats time durations given in milliseconds."""
    int_time = int(a)
    
    days =  '' if a < 86400 else f"{int_time // 86400}d "
    hours = '' if a < 3600 else f"{(int_time % 86400) // 3600}h "
    mins = '' if a < 60 else f"{(int_time % 3600) // 60}m "
    secs = f"{a % 60:.2f}s"

    return f"{days}{hours}{mins}{secs}"

class ProcessTimer:
    """A class for measuring process runtimes.
    
    Use the `start()` and `stop()` methods.

    Args:
        verbose (bool, optional): If True (default), prints the process time
            when stopped.
    """

    def __init__(self, verbose=True):
        self.proc_name = None
        self.verbose = verbose
        self.start_time = None
        self.add_time = 0

    def start(self, proc_name=None):
        """Starts the timer. If the timer is already running, stops it first."""

        if self.start_time: ### If the timer is already running
            self.stop()

        self.start_time = time.time()
        self.proc_name = proc_name

    def pause(self):
        """Pauses the timer. Can be continued with `start()`"""

        if self.start_time: ### Only if the timer is already running
            self.add_time += time.time()-self.start_time
            self.start_time = None

    def stop(self, proc_name=None):
        """Stops and clears the timer and returns the process name and run time."""
        
        ### Only if the timer is already running or was paused
        if self.start_time or self.add_time:

            if proc_name:
                self.proc_name = proc_name

            proc_time = (time.time()-self.start_time if self.start_time else 0) + self.add_time
            return_value = f"{self.proc_name+': ' if self.proc_name else ''}{format_time(proc_time)}"
            self.start_time = None
            self.proc_name = None

            if self.verbose:
                print(return_value)

            return return_value

def check_path(path, mode=0):
    """Checks if path exists. If it doesn't, reacts based on mode.
    
    There are three modes available:
        0 - Raise error if path does not exist.
        1 - Create path (if directory) and raise error if parents are missing.
        2 - Create path (if directory) with all its missing parents.

    Args:
        path (path): The path to check.
        mode (int): Decides what happens if the path is missing. Can be 0, 1, 2.

    Returns:
        None
    """

    path = Path(path)

    if mode == 0:
        if not path.exists():
            raise FileNotFoundError(f"This path does not exist: {path}")
    
    elif mode == 1:
        if not path.exists():
            path.mkdir()

    elif mode == 2:
        if not path.exists():
            path.mkdir(parents=True)

    else:
        raise ValueError("Mode has to be 0, 1 or 2")

def load_files(paths, pattern='*.*', recursive=False, str_out=False):
    """Creates a list of input files.

    The paths can be:
        - Path to a single file matching the pattern.
        - Path to a directory containing matching files.
        - Python list of paths to files and directories.
        - Path to a .txt file containing a list with each entry in a new line.

    Args:
        paths (path or list): Input path or list.
        pattern (str, optional): The glob-style pattern the paths must match.
        recursive (bool, optional): Whether to scan directories recursively.
        str_out (bool, optional): If True, returns strings instead of pathlib
            Path objects.

    Returns:
        file_list (list): A list of paths. Only includes existing paths.
    """

    prefix = '**/' if recursive else '' ### Decides if globbing is recursive or not

    ### Creating a list if a simple path is given
    if isinstance(paths, (str, Path)):
        paths = Path(paths)

        ### If it is a txt, treat it as a list
        if paths.match('*.txt'):
            with open(paths, 'rt') as f:
                paths = [Path(i) for i in [line.split('#', 1)[0].strip().strip('\u0022\u0027') for line in f] if i]
                ### Empty lines are omitted, anything after a hash (#) is discarded.

        ### If it is another file format or a directory, wrap it into a list
        else:
            paths = [paths]

    ### Processing the list
    if isinstance(paths, list):
        file_list = []

        for i in paths:
            i = Path(i)
            if i.exists(): ### Nonexistent paths are simply not included

                ### Directories are globbed
                if i.is_dir():
                    file_list += list(i.glob(f"{prefix}{pattern}"))
                
                ### Matching paths are added to the list
                elif i.match(pattern):
                    file_list += [i]
        
        return file_list if not str_out else [str(i) for i in file_list]

    else:
        raise TypeError("Paths need to match the pattern or be a list or directory.")
    
### File Container

class DirInterface:
    """File container interface for accessing directories."""
    def __init__(self, cont_path):
        self.container = Path(cont_path)
    
    def filter(self, pattern):
        """Returns a list of files that fit the pattern."""
        return list(self.container.rglob(pattern))

    def get(self, filename):
        """Returns a path or file-like object to work on."""
        if Path(filename).exists():
            return filename
        else:
            raise FileNotFoundError(f"{filename} does not exist.")

    def extract(self, filename, out_path):
        """Copies a file to another location."""
        shutil.copy(filename, out_path)

    def close(self):
        pass

class ZIPInterface:
    """File container interface for accessing .zip files."""
    def __init__(self, cont_path):
        self.container = zipfile.ZipFile(cont_path, 'r')
    
    def filter(self, pattern):
        """Returns a list of files that fit the pattern."""
        return fnmatch.filter(self.container.namelist(), pattern)

    def get(self, filename):
        """Returns a path or file-like object to work on."""
        return self.container.open(filename)

    def extract(self, filename, out_path):
        """Copies a file to another location."""
        zip_info = self.container.getinfo(filename)
        out_path = Path(out_path)
        zip_info.filename = out_path.name
        self.container.extract(zip_info, out_path.parent)

    def close(self):
        self.container.close()

class TARInterface:
    """File container interface for accessing .tar files."""
    def __init__(self, cont_path):
        self.container = tarfile.open(cont_path, 'r')
    
    def filter(self, pattern):
        """Returns a list of files that fit the pattern."""
        return fnmatch.filter(self.container.getnames(), pattern)

    def get(self, filename):
        """Returns a path or file-like object to work on."""
        return self.container.extractfile(filename)

    def extract(self, filename, out_path):
        """Copies a file to another location."""
        tar_info = self.container.getmember(filename)
        out_path = Path(out_path)
        tar_info.name = out_path.name
        self.container.extract(tar_info, out_path.parent)

    def close(self):
        self.container.close()
    
class FileContainer:
    """Context manager to provide standardize handling of file containers.

    Use the `with` statement.
    
    Args:
        cont_path (path): Path to the file container.
        cont_type (str, optional): Type of the container (dir, zip or tar).
            Attempts to detect it if not provided.
    """
    
    def __init__(self, cont_path, cont_type=None):

        cont_path = Path(cont_path)
        
        if cont_type:
            subclasses = {
                'zip': ZIPInterface,
                'tar': TARInterface,
                'dir': DirInterface,
            }

            if cont_type not in subclasses.keys():
                raise ValueError("Unknown container type.")

            self.interface = subclasses[cont_type](cont_path)
        
        else:
            if cont_path.is_dir():
                self.interface = DirInterface(cont_path)

            else:
                suffix = cont_path.suffix

                if suffix == '.zip':
                    self.interface = ZIPInterface(cont_path)

                elif suffix in ['.tar', '.tgz', '.tbz', '.txz', '.gz', '.bz2', '.xz']:
                    self.interface = TARInterface(cont_path)

                else:
                    raise ValueError("Unknown container type.")

    def __enter__(self):
        return self.interface

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.interface.close()
        if exc_type is not None:
            return False