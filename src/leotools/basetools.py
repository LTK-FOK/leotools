""""This module contains tools supplmentary to geoinformatic applications."""

from pathlib import Path
import time
import fnmatch

### Basic stuff

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
            int_proc_time = int(proc_time)
            
            days =  '' if proc_time < 86400 else f"{int_proc_time // 86400}d "
            hours = '' if proc_time < 3600 else f"{(int_proc_time % 86400) // 3600}h "
            mins = '' if proc_time < 60 else f"{(int_proc_time % 3600) // 60}m "
            secs = f"{proc_time % 60:.2f}s"

            return_value = f"{self.proc_name+': ' if self.proc_name else ''}{days}{hours}{mins}{secs}"
            self.start_time = None
            self.proc_name = None

            if self.verbose:
                print(return_value)

            return return_value

def check_path(path, mode=0):
    """Checks if path exists. If it doesn't, reacts based on mode.
    
    There are three modes available:
        0 - Raise error if path does not exist.
        1 - Create path and raise error if parents are missing.
        2 - Create path with all its missing parents.

    Args:
        path (path or str): The path to check.
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
            if i.exists(): ### Nonexistent paths are simply not included
                i = Path(i)

                ### Directories are globbed
                if i.is_dir():
                    file_list += list(i.glob(f"{prefix}{pattern}"))
                
                ### Matching paths are added to the list
                elif i.match(pattern):
                    file_list += [i]
        
        return file_list if not str_out else [str(i) for i in file_list]

    else:
        raise TypeError("Paths need to match the pattern or be a list or directory.")