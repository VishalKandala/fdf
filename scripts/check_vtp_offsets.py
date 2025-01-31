#!/usr/bin/env python3

import re
import struct

filename = "velocity00000.vtp"

with open(filename, 'rb') as f:
    content = f.read()

# 1) Split file content at <AppendedData encoding="raw"> and </AppendedData>
try:
    header_part, appended_part = content.split(b'<AppendedData encoding="raw">', 1)
    appended_data, _ = appended_part.split(b'</AppendedData>', 1)
except ValueError:
    print("Could not find <AppendedData encoding=\"raw\"> or </AppendedData> in file.")
    exit(1)

# 2) Extract declared offsets from the XML header
offsets = re.findall(rb'offset="(\d+)"', header_part)
offsets = [int(o) for o in offsets]
print("Found offsets in header:", offsets)

# 3) Step through each declared offset, read its size block, and verify
#    By convention, appended data often starts with a single '_' (0x5F) before the first 4-byte size.
pos = 0
if appended_data[:1] == b'_':  
    pos = 1  # skip the leading underscore

for i, declared_offset in enumerate(offsets):
    # Check if our current 'pos' matches the declared offset
    if pos != declared_offset:
        print(f"WARNING: Declared offset={declared_offset} but actual pos={pos} for DataArray index={i}")

    # Read the 4-byte block size (little-endian uint32)
    block_size = struct.unpack('<I', appended_data[pos:pos+4])[0]
    pos += 4

    # Move past the actual data payload
    pos += block_size

    print(f"DataArray {i}: read block_size={block_size}, next position={pos}")

print("Finished checking appended data vs. declared offsets.")
