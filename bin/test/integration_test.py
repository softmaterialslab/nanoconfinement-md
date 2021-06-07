import pytest
import subprocess
import sys
import threading
import pandas as pd
import math

def exec_command(command, cwd=None):
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               universal_newlines=True)
    while True:
        output = process.stdout.readline()
        print(output.strip())
        # Do something else
        return_code = process.poll()
        if return_code is not None:
            print('RETURN CODE', return_code)
            # Process has finished, read rest of the output
            for output in process.stdout.readlines():
                print(output.strip())
            break


def test_example_one():
    exec_command(['md_simulation_confined_ions', '-T', '0.001', '-Z', '3', '-p', '1', '-n', '-1', '-c', '.5', '-d', '0.474', '-a', '0.627', '-i', '0', '-S', '10000'], cwd="..")
    df = pd.read_csv('outfiles/energy.dat', sep=r"\s+", header=None)
    t = df.iat[-1, 6]
    assert math.isclose(t, 2310.45, rel_tol=1e-1)






