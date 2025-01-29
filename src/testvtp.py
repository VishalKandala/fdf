import struct
import xml.etree.ElementTree as ET
import sys
import os

def validate_and_fix_vtp(filename):
    print(f"\n🔍 Validating and fixing VTP File: {filename}\n")

    if not os.path.exists(filename):
        print(f"❌ ERROR: File {filename} not found!")
        return

    try:
        with open(filename, "rb") as f:
            data = f.read()

        # Extract XML header
        header_end = data.find(b'<AppendedData encoding="raw">')
        if header_end == -1:
            print("❌ ERROR: Missing or invalid <AppendedData> section.")
            return

        xml_header = data[:header_end].decode("utf-8", errors="replace") + "</AppendedData></VTKFile>"

        # Validate XML structure
        try:
            root = ET.fromstring(xml_header)
            print("✅ XML header is valid.")
        except ET.ParseError as e:
            print(f"❌ XML parsing error: {e}. Attempting to fix...")
            xml_header = xml_header.strip()
            if not xml_header.endswith("</VTKFile>"):
                xml_header += "</AppendedData></VTKFile>"
            try:
                root = ET.fromstring(xml_header)
                print("✅ XML header fixed and valid.")
            except ET.ParseError as e:
                print(f"❌ Critical XML parsing error: {e}. Cannot proceed.")
                return

        # Locate binary section
        binary_start = data.find(b'<AppendedData encoding="raw">') + len("<AppendedData encoding=\"raw\">\n_")
        binary_data = data[binary_start:]

        if not binary_data:
            print("❌ ERROR: Missing binary data section.")
            return

        # Validate and fix binary blocks
        def read_block(data, offset, name, dtype, size_per_element):
            if offset >= len(data):
                print(f"❌ ERROR: {name} block offset out of bounds.")
                return None, offset

            try:
                block_size = struct.unpack("<I", data[offset:offset + 4])[0]
                offset += 4
                block_data = data[offset:offset + block_size]
                offset += block_size

                if len(block_data) != block_size:
                    print(f"❌ ERROR: {name} block size mismatch. Expected {block_size}, got {len(block_data)}.")
                    return None, offset

                num_elements = block_size // size_per_element
                unpacked_data = struct.unpack(f"<{num_elements}{dtype}", block_data)
                print(f"✅ {name} block validated: {block_size} bytes, {num_elements} elements.")
                return unpacked_data, offset
            except (struct.error, ValueError) as e:
                print(f"❌ ERROR: Failed to read {name} block: {e}")
                return None, offset

        offset = 0
        position_data, offset = read_block(binary_data, offset, "Position", "d", 8)
        velocity_data, offset = read_block(binary_data, offset, "Velocity", "d", 8)
        connectivity_data, offset = read_block(binary_data, offset, "Connectivity", "i", 4)
        offsets_data, offset = read_block(binary_data, offset, "Offsets", "i", 4)

        # Rebuild the fixed file
        fixed_filename = filename.replace(".vtp", "_fixed.vtp")
        with open(fixed_filename, "wb") as f:
            # Write fixed XML header
            f.write(xml_header.encode("utf-8"))
            f.write(b"\n<AppendedData encoding=\"raw\">\n_")

            # Write fixed binary blocks
            def write_block(f, name, data, dtype, size_per_element):
                if data is None:
                    print(f"⚠️ WARNING: Skipping {name} block due to errors.")
                    return
                block_size = len(data) * size_per_element
                f.write(struct.pack("<I", block_size))
                f.write(struct.pack(f"<{len(data)}{dtype}", *data))
                print(f"✅ {name} block written: {block_size} bytes.")

            write_block(f, "Position", position_data, "d", 8)
            write_block(f, "Velocity", velocity_data, "d", 8)
            write_block(f, "Connectivity", connectivity_data, "i", 4)
            write_block(f, "Offsets", offsets_data, "i", 4)

            # Close AppendedData
            f.write(b"\n</AppendedData>\n</VTKFile>")

        print(f"✅ Fixed VTP file saved as: {fixed_filename}")

    except Exception as e:
        print(f"❌ Critical error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 vtp_build.py <filename>")
    else:
        validate_and_fix_vtp(sys.argv[1])
