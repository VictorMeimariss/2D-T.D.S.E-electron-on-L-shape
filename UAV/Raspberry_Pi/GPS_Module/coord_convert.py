# Class for converting from GPS Module format to Ardupilot GPAA format
class Convert:
    def __init__(self, buff):
        self.buff = buff

    def parse_cgnssinfo(self, buff):
        # Strip leading/trailing whitespace and hidden characters
        buff = buff.strip()

        if not buff.startswith('+CGNSSINFO:'):
            return None

        parts = buff.split(": ")[1].split(",")

        if len(parts) < 16:
            return None

        fix_status = parts[3]
        used_sv = parts[2]
        utc_time = parts[9]
        latitude = parts[4]
        lat_dir = parts[5]
        longitude = parts[6]
        lon_dir = parts[7]
        altitude = parts[10]
        hdop = parts[13]

        quality = "1" if fix_status in ["01", "02"] else "0"
        gga = f"$GPGGA,{utc_time},{latitude},{lat_dir},{longitude},{lon_dir},{quality},{used_sv},{hdop},{altitude},M,0.0,M,,"
        checksum = self.checksum(gga)
        gga += f"*{checksum}"
        return gga

    def checksum(self, sentence):
        checksum_value = 0
        for char in sentence[1:]:
            if char == '*':
                break
            checksum_value ^= ord(char)
        return f"{checksum_value:02X}"