import platform

class Settings:
    if platform.system() == "Linux":
        SERIAL_DEVICE = "/dev/ttyUSB0"
    elif platform.system() == "Darwin":
        SERIAL_DEVICE = "/dev/tty.usbserial-0001"
    else:
        raise Exception("OS not supported")

