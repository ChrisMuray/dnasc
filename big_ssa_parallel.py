#!/usr/bin/python3

from concurrent.futures import ProcessPoolExecutor
import subprocess

def sim():
    subprocess.call(['./con.py'])

if __name__ == '__main__':
    nsims = 10
    with ProcessPoolExecutor(max_workers=nsims) as executor:
        subprocess_futures = [executor.submit(sim) for i in range(nsims)]

