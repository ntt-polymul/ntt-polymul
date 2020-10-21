#!/usr/bin/env python3
import subprocess
import sys
import serial
import hashlib
from config import Settings

def test(scheme, refHash):
    subprocess.check_call("make clean", shell=True)
    binary = f"bin/crypto_kem_{scheme}_m4_testvectors.bin"
    make = f"make IMPLEMENTATION_PATH=crypto_kem/{scheme}/m4 {binary}"
    subprocess.check_call(make, shell=True)

    try:
        subprocess.check_call(f"st-flash write {binary} 0x8000000", shell=True)
    except:
        print("flashing failed --> retry")
        return test(scheme, refHash)

    # get serial output and wait for '#'
    with serial.Serial(Settings.SERIAL_DEVICE, 115200, timeout=10) as dev:
        log = b""
        while True:
            device_output = dev.read()
            if device_output == b'':
                print("timeout --> retry")
                return test(scheme, refHash)
            sys.stdout.buffer.write(device_output)
            sys.stdout.flush()
            log += device_output
            if device_output == b'#':
                break

    log = log.decode(errors="ignore")
    log = log.split("=")[-1].split("#")[0][1:]
    hash = hashlib.sha3_256(log.encode("utf-8")).hexdigest()
    assert hash == refHash.lower()

test("lightsaber", "34D3126F8AAEF14DAB8413F618E2734C04DD74EC5CA1C94D2A4A0E2D3EA718F2")
test("saber",      "8CBA64557814372316871038280A5195A62511EEBC4792F802AE425F82158F96")
test("firesaber",  "C8955448B1F7AABBDFE80ADB924201154D324A8790E26E3B15BACDBEED8CF291")

test("ntruhrss701","3E544517228C8CD0B15A37B2DC573C46F133B552FEB1B69329A15BC2C8327C05")
test("ntruhps2048509","ED02BDBA3669EDE66F98067A056DD5355095215F0A2EF01CBB59C1C77CDB9447")
test("ntruhps2048677","543293E03C0BE5DF3E5FDC745008041075696A8825D12D21CBFBC8B86128952D")
test("ntruhps4096821","F9687F87550259DFDA3838699CD81F756D04E899E346051B75AD6C2DFA4289CE")

print("all tests passed.")
