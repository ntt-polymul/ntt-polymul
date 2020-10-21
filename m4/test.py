#!/usr/bin/env python3
import subprocess
import sys
import serial
from config import Settings


def test(scheme):
    subprocess.check_call("make clean", shell=True)
    binary = f"bin/crypto_kem_{scheme}_m4_test.bin"
    make = f"make IMPLEMENTATION_PATH=crypto_kem/{scheme}/m4 {binary}"
    subprocess.check_call(make, shell=True)

    try:
        subprocess.check_call(f"st-flash write {binary} 0x8000000", shell=True)
    except:
        print("flashing failed --> retry")
        return test(scheme)
    # get serial output and wait for '#'
    with serial.Serial(Settings.SERIAL_DEVICE, 115200, timeout=10) as dev:
        log = b""
        while True:
            device_output = dev.read()
            if device_output == b'':
                print("timeout --> retry")
                return test(scheme)
            sys.stdout.buffer.write(device_output)
            sys.stdout.flush()
            log += device_output
            if device_output == b'#':
                break

    log = log.decode(errors="ignore")
    assert log.count("ERROR") == 0 and log.count("OK") == 30

test("lightsaber")
test("saber")
test("firesaber")

test("ntruhrss701")
test("ntruhps2048509")
test("ntruhps2048677")
test("ntruhps4096821")
print("all tests passed.")
