#!/usr/bin/env python3
import sys
import re
import xml.etree.ElementTree as ET
import struct
import os

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <file.vtp>")
        sys.exit(1)

    filename = sys.argv[1]
    if not os.path.isfile(filename):
        print(f"[ERROR] File '{filename}' does not exist.")
        sys.exit(1)

    with open(filename, 'rb') as f:
        file_bytes = f.read()

    ############################################
    # 1) Parse XML structure (up to AppendedData)
    ############################################
    # Try to split off the <AppendedData> section in raw form
    try:
        # We expect something like '<AppendedData encoding="raw">_... </AppendedData>'
        content_str = file_bytes.decode('utf-8', errors='replace')  # decode for XML
    except UnicodeDecodeError:
        print("[ERROR] Failed to decode file as UTF-8 to parse XML. Possibly corrupted.")
        sys.exit(1)

    # We'll search for <AppendedData encoding="raw"> and </AppendedData>
    # Then parse everything before <AppendedData ...> as XML.
    appended_start = content_str.find('<AppendedData encoding="raw">')
    appended_end   = content_str.find('</AppendedData>')
    if appended_start == -1 or appended_end == -1:
        print("[ERROR] Could not locate <AppendedData ...> or </AppendedData> tags in file.")
        sys.exit(1)

    header_part = content_str[:appended_start]
    appended_part = content_str[appended_start:appended_end]

    # Basic XML parse on header_part plus a dummy close tag
    # because we truncated the real <AppendedData> element.
    # We'll artificially patch it so ElementTree can parse it.
    # Re-inject a minimal well-formed close tag.
    testable_xml = header_part + "<AppendedData></AppendedData></VTKFile>"
    try:
        root = ET.fromstring(testable_xml)
    except ET.ParseError as e:
        print(f"[ERROR] XML parsing failed for the header portion: {e}")
        sys.exit(1)

    # Check top-level <VTKFile> attributes
    if root.tag != 'VTKFile':
        print(f"[ERROR] Top-level element is '{root.tag}', expected 'VTKFile'.")
    file_type = root.attrib.get('type', '')
    if file_type != 'PolyData':
        print(f"[WARN] <VTKFile type='{file_type}'> but this script expects 'PolyData'.")
    byte_order = root.attrib.get('byte_order', '')
    print(f"[INFO] VTKFile type='{file_type}', byte_order='{byte_order}'")

    # Go inside <PolyData> -> <Piece>
    polydata_elem = root.find('PolyData')
    if polydata_elem is None:
        print("[ERROR] No <PolyData> element found.")
        sys.exit(1)
    piece_elem = polydata_elem.find('Piece')
    if piece_elem is None:
        print("[ERROR] No <Piece> element found inside <PolyData>.")
        sys.exit(1)

    # Gather <Piece> attributes
    npoints_str = piece_elem.attrib.get('NumberOfPoints', '0')
    nverts_str  = piece_elem.attrib.get('NumberOfVerts', '0')
    try:
        npoints = int(npoints_str)
        nverts  = int(nverts_str)
    except ValueError:
        print("[ERROR] <Piece> NumberOfPoints/NumberOfVerts not valid integers.")
        sys.exit(1)
    print(f"[INFO] Piece has NumberOfPoints={npoints}, NumberOfVerts={nverts}.")

    # Collect all <DataArray> with 'format="appended"' and offsets
    data_arrays = piece_elem.findall('.//DataArray[@format="appended"]')
    if not data_arrays:
        print("[ERROR] No <DataArray format=\"appended\"> elements found; cannot check offsets.")
        sys.exit(1)

    # Each found DataArray => read 'offset' attribute
    # Keep them in the order encountered
    array_info = []
    for da in data_arrays:
        offset_str = da.attrib.get('offset', None)
        if offset_str is None:
            print("[ERROR] <DataArray> missing 'offset'.")
            continue
        try:
            offset_val = int(offset_str)
        except ValueError:
            print(f"[ERROR] <DataArray offset=\"{offset_str}\"> is not an integer.")
            continue

        name = da.attrib.get('Name', '')
        dtype = da.attrib.get('type', '')
        comps_str = da.attrib.get('NumberOfComponents', '1')
        try:
            comps = int(comps_str)
        except ValueError:
            comps = 1

        array_info.append({
            'name': name,
            'offset': offset_val,
            'type': dtype,
            'components': comps
        })

    # Print summary of DataArrays
    print(f"[INFO] Found {len(array_info)} appended DataArray(s).")
    for i, info in enumerate(array_info):
        print(f"  - DataArray {i} Name='{info['name']}', offset={info['offset']}, "
              f"type='{info['type']}', comps={info['components']}")

    #######################################
    # 2) Inspect the raw appended data area
    #######################################
    # In the actual byte content, find the start of "<AppendedData encoding=\"raw\">"
    # and skip past it. Typically there's an underscore '_' then raw data.
    # We'll find exactly the substring "<AppendedData encoding=\"raw\">" and measure from there.

    appended_data_str = "<AppendedData encoding=\"raw\">"
    start_pos = file_bytes.find(appended_data_str.encode('utf-8'))
    if start_pos < 0:
        print("[ERROR] Could not locate raw appended data region in bytes.")
        sys.exit(1)

    # Move start_pos to the end of that tag
    start_pos += len(appended_data_str)

    # We expect the first character to be '_'
    # Some code might include a newline after '_'; let's detect both.
    # We'll just skip any whitespace between '>' and '_'.
    # Then appended binary starts exactly at offset 0 relative to that next byte.
    raw_section = file_bytes[start_pos:]

    # Trim leading whitespace
    raw_section = raw_section.lstrip()
    if not raw_section:
        print("[ERROR] No bytes found after <AppendedData encoding=\"raw\">; file truncated?")
        sys.exit(1)

    if raw_section[0:1] != b'_':
        print("[WARN] The raw section does not begin with underscore '_'. "
              "This can break appended format. Found:", raw_section[0:1])
    else:
        # Chop off the underscore itself
        raw_section = raw_section[1:]
        # Now offset=0 means the first byte of raw_section.

    # We'll also look for any sign of newline still stuck after the underscore
    if len(raw_section) > 0 and raw_section[0] in (10, 13):
        print("[WARN] There's a newline or CR right after '_'. "
              "That can cause offset misalignment for some VTK readers.")

    # For thoroughness, let's see if the last chunk includes the closing </AppendedData>
    # We'll rely on appended_end to confirm we found it earlier. Not re-checking here.

    ############################################
    # 3) Compare declared offsets vs actual data
    ############################################
    # Typically, each appended block is:
    #   4-byte length (little-endian), then 'length' bytes of data.
    # We'll see if that pattern holds. We'll do:
    #
    #   position=offset_0
    #   read 4 bytes => block_size => skip block_size => next offset
    #
    #   compare next offset with offset_1 ...
    #
    # This is to confirm each offset is correct. If there's no 4-byte block size
    # actually present, we might see nonsense or run out of data.

    pos = 0
    for i, info in enumerate(array_info):
        name = info['name']
        declared_offset = info['offset']

        # Check if our 'pos' matches the declared offset
        if pos != declared_offset:
            print(f"[WARN] DataArray '{name}' mismatch: declared offset={declared_offset}, "
                  f"but actual pos={pos}. (Possible leftover newlines or missing block-size headers?)")

        # Try to read the 4-byte block size
        if pos + 4 > len(raw_section):
            print(f"[ERROR] Not enough bytes to read the 4-byte size header for DataArray '{name}'.")
            break

        block_size_bytes = raw_section[pos:pos+4]
        block_size = struct.unpack('<I', block_size_bytes)[0]
        pos += 4  # move past the size field

        # Check that block_size won't exceed total appended data
        end_of_block = pos + block_size
        if end_of_block > len(raw_section):
            print(f"[ERROR] DataArray '{name}' declared block_size={block_size} "
                  f"overruns appended data boundary.")
            break

        # Move pos to the end of this block
        pos = end_of_block
        print(f"[INFO] DataArray '{name}': size={block_size}, next pos={pos}")

    # We'll do a final check if we ended near the end
    leftover = len(raw_section) - pos
    if leftover > 0:
        print(f"[INFO] {leftover} leftover bytes remain after the last declared DataArray. "
              "This can be normal or might indicate additional data not referenced by <DataArray> offsets.")

    leftover_bytes = raw_section[pos:]
    print("[INFO] Leftover bytes (hex):", leftover_bytes.hex(' ', -4))
    print("[INFO] Leftover bytes (ASCII):", leftover_bytes.decode('ascii', errors='replace'))


    print("[DONE] Finished checks.")


if __name__ == "__main__":
    main()
