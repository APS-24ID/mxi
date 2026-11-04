"""A command run in the background, its output a line at a time: no tkinter
here -- the window polls lines() from its own loop."""

import os
import queue
import shlex
import signal
import subprocess
import threading


def as_text(argv):
    """The command as one would type it, for showing and copying."""
    return shlex.join(argv)


class Runner:
    """Runs one command at a time. lines() returns what it has printed since
    the last call; finished() is its exit code once it has ended, else None."""

    def __init__(self):
        self._process = None
        self._lines = queue.Queue()
        self._code = None
        self._thread = None

    def start(self, argv, cwd=None, env=None):
        if self.running():
            raise RuntimeError("a command is running already")
        self._code = None
        self._process = subprocess.Popen(
            argv,
            cwd=cwd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            text=True,
            bufsize=1,
            # Its own process group, so that cancel stops whatever it started.
            start_new_session=(os.name == "posix"),
        )
        self._thread = threading.Thread(target=self._read, daemon=True)
        self._thread.start()

    def _read(self):
        process = self._process
        for line in process.stdout:
            self._lines.put(line.rstrip("\n").split("\r")[-1])
        process.stdout.close()
        self._code = process.wait()

    def running(self):
        return self._process is not None and self._code is None

    def lines(self):
        out = []
        while True:
            try:
                out.append(self._lines.get_nowait())
            except queue.Empty:
                return out

    def finished(self):
        """The exit code once the command and its output are done, else None."""
        if self._thread is not None and self._thread.is_alive():
            return None
        return self._code

    def cancel(self):
        if not self.running():
            return
        try:
            if os.name == "posix":
                os.killpg(self._process.pid, signal.SIGTERM)
            else:
                self._process.terminate()
        except ProcessLookupError:
            pass
