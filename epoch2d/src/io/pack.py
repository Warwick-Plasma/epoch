#! /usr/bin/env python

import binascii
import struct
import os
import sys
import time
import subprocess as sp
import tarfile
import hashlib
import platform
import gzip

prefix="sdf"
pack_source_code = True
pack_git_diff = True
pack_git_diff_from_origin = True
generate_checksum = True

archive="source_info_archive.tgz"
hexdump="source_info_hexdump.txt"
gitdiff="source_info_gitdiff.txt"
varname="%s_bytes" % prefix
diffname="%s_diff_bytes" % prefix
module_name="%s_source_info" % prefix
outfile=sys.argv[1]
incfile=outfile.split('.')[0] + '_include.inc'

f77_output = False
nbytes = 8
nelements = 0
padding = 0
vname=varname

if f77_output:
  inc_handle=open(incfile,"w")
  linestart=6*' '
  linecont=5*' '+'&'
  suffix=''
  ncolumns=72
  ncontinuation=19
else:
  linestart=''
  linecont=''
  suffix='&'
  ncolumns=132 # gfortran ignores the F90 standard of 139
  ncontinuation=39


def byteswap4(s):
  s = binascii.unhexlify(s)
  a, = struct.unpack('>L',s)
  s = struct.pack('<L',a)
  return binascii.hexlify(s)


def byteswap8(s):
  s = binascii.unhexlify(s)
  a, = struct.unpack('>Q',s)
  s = struct.pack('<Q',a)
  return binascii.hexlify(s)


def wrapped(string):
  global of, linestart, ncolumns
  ostring = linestart + string
  rem = len(ostring)
  if f77_output:
    while (rem > ncolumns):
      of.write(ostring[:ncolumns]+'\n')
      ostring = linecont + ostring[ncolumns:]
      rem = len(ostring)
  else:
    while (rem > ncolumns):
      of.write(ostring[:ncolumns-1]+'&\n')
      ostring = '&' + ostring[ncolumns-1:]
      rem = len(ostring)
  if rem > 0:
    of.write(ostring)
  of.write('\n')


def print_character(name,value):
  global vname,of
  #ilen = len(value)
  ilen = 256
  var="%s_%s" % (vname,name)
  if f77_output:
    of=inc_handle
    wrapped("CHARACTER*%i %s" % (ilen,var))
    wrapped("COMMON/c_%s/%s" % (vname,var))
    of=out_handle
    wrapped("CHARACTER*%i %s" % (ilen,var))
    wrapped("COMMON/c_%s/%s" % (vname,var))
    wrapped("DATA %s/'%s'/" % (var,value))
  else:
    ilen = len(value)
    if ilen == 0: ilen = 1
    wrapped("CHARACTER(LEN=%i) :: %s = '%s'" % (ilen,var,value))


def print_integer(name,value):
  global vname,of
  var="%s_%s" % (vname,name)
  if f77_output:
    of=inc_handle
    wrapped("INTEGER " + var)
    wrapped("COMMON/i_%s/%s" % (vname,var))
    of=out_handle
    wrapped("INTEGER " + var)
    wrapped("COMMON/i_%s/%s" % (vname,var))
    wrapped("DATA %s/%i/" % (var,value))
  else:
    wrapped("INTEGER, PARAMETER :: %s = %i" % (var,value))


def get_bytes_checksum(files):
  global checksum_type
  if not generate_checksum:
    checksum_type = ''
    return ''
  cksum = hashlib.new('sha256')
  for name in files:
    f = open(name)
    while True:
      data = f.read(cksum.block_size)
      if not data:
        break
      cksum.update(data)
  checksum_type = 'sha256'
  return cksum.hexdigest()


def write_data_bytes(filename, varname):
  global mimetype, of
  global linestart, linecont, suffix, ncolumns, ncontinuation

  f=open(filename)
  d=f.read()
  dhex=d.encode('hex')
  f.close()
  os.remove(filename)

  nelements = (len(d)+nbytes-1) / nbytes
  padding = nelements * nbytes - len(d)
  dhex += '00' * padding

  vname=varname
  print_character('mimetype', mimetype)
  print_integer('padding', padding)
  print_integer('len', nelements)
  print_integer_array(nelements)

  nwidth = len("z'',") + 2 * nbytes
  nper_line_body = (ncolumns - 1) / nwidth
  sdata = linestart + "DATA(%s(i),i=%i,%i)/" % (varname,nelements,nelements)
  nper_line_first = (ncolumns - len(sdata) - 1) / nwidth
  nper_segment = nper_line_first + nper_line_body * ncontinuation

  i0 = 0
  segline = 0
  elements_written = 0
  while elements_written < nelements:
    ss = ""
    if segline == 0:
      i1 = min(i0 + nper_segment, nelements)
      ss += linestart + "DATA(%s(i),i=%i,%i)/" % (varname,i0+1,i1)
      i0 = i1
    else:
      ss += linecont

    shex = dhex[2*nbytes*elements_written:2*nbytes*(elements_written+1)]
    if nbytes == 4:
      shex = byteswap4(shex)
    elif nbytes == 8:
      shex = byteswap8(shex)
    ss += "z'%s'" % shex
    elements_written = elements_written + 1
    if segline == 0:
      nper_line = nper_line_first - 1
    else:
      nper_line = nper_line_body - 1
    n = 0
    while n < nper_line and elements_written != nelements:
      shex = dhex[2*nbytes*elements_written:2*nbytes*(elements_written+1)]
      if nbytes == 4:
        shex = byteswap4(shex)
      elif nbytes == 8:
        shex = byteswap8(shex)
      ss += ",z'%s'" % shex
      elements_written = elements_written + 1
      n = n + 1

    if elements_written == nelements or segline == ncontinuation:
      ss += "/\n"
    else:
      ss += "," + suffix + "\n"

    of.write(ss)

    if segline == ncontinuation:
      segline = 0
    else:
      segline = segline + 1


def print_integer_array(value):
  global of,nbytes,vname
  if value == 0: value = 1
  if f77_output:
    of=inc_handle
    wrapped("INTEGER*%i %s(%s_len)" % (nbytes,vname,vname))
    of=out_handle
    wrapped("INTEGER*%i %s(%i)" % (nbytes,vname,value))
  else:
    wrapped("INTEGER(%i) :: %s(%i)" % (nbytes,vname,value))


fnull=open(os.devnull,'w')
try:
  git_version = sp.check_output("git describe --always --long --dirty",
                  shell=True,stderr=fnull).rstrip()
except:
  git_version = ''
  pack_git_diff = False

tsec = time.time()
compile_date = int(round(tsec))
compile_date_string = time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime(tsec))
compile_machine_info = ' '.join((platform.node(),platform.platform()))
compiler_info=sys.argv[2]
compiler_flags=sys.argv[3]

#filelist = sp.check_output("git ls-files --cached --no-empty-directory "
#                + "--full-name", shell=True, stderr=fnull).rstrip()

filelist = sys.argv[4:]
if filelist == []:
  pack_source_code = False


# Now write the file

out_handle=open(outfile,"w")
of=out_handle

if f77_output:
  of.write(linestart + "SUBROUTINE %s_source_info\n" % module_name)
else:
  of.write(linestart + "MODULE %s\n\n" % module_name)
  of.write(linestart + "IMPLICIT NONE\n\n")

print_character('git_version', git_version)
print_character('compile_date_string', compile_date_string)
print_character('compile_machine_info', compile_machine_info)
print_character('compiler_info', compiler_info)
print_character('compiler_flags', compiler_flags)
print_integer('compile_date', compile_date)

if pack_source_code or pack_git_diff:
  if f77_output:
    of.write(linestart + "INTEGER i\n\n")
  else:
    of.write(linestart + "INTEGER, PRIVATE :: i\n")

vname = varname
checksum_type = ''
checksum = ''
if filelist != []:
  checksum = get_bytes_checksum(filelist)
print_character('checksum_type', checksum_type)
print_character('checksum', checksum)
if not pack_source_code:
  mimetype = ''
  print_character('mimetype', mimetype)
  print_integer('padding', padding)
  print_integer('len', 0)
  print_integer_array(0)
else:
  tar = tarfile.open(archive, "w:gz")
  for name in filelist:
    tar.add(name)
  tar.close()
  mimetype = 'application/x-tar-gz'

  write_data_bytes(archive, vname)


vname = diffname
checksum_type = ''
checksum = ''
if not pack_git_diff:
  mimetype = ''
  print_character('checksum_type', checksum_type)
  print_character('checksum', checksum)
  print_character('mimetype', mimetype)
  print_integer('padding', padding)
  print_integer('len', 0)
  print_integer_array(0)
else:
  if pack_git_diff_from_origin:
    sp.call(["git diff origin/master > %s" % gitdiff], shell=True)
  else:
    sp.call(["git diff > %s" % gitdiff], shell=True)
  if os.path.getsize(gitdiff) != 0:
    checksum = get_bytes_checksum([gitdiff])

    zgitdiff = gitdiff + '.gz'
    f_in  = open(gitdiff, 'rb')
    f_out = gzip.open(zgitdiff, 'wb')
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(gitdiff)
    os.rename(zgitdiff,gitdiff)
  mimetype = 'application/x-gzip'

  print_character('checksum_type', checksum_type)
  print_character('checksum', checksum)

  write_data_bytes(gitdiff, vname)

if f77_output:
  of.write(linestart + "END SUBROUTINE\n")
else:
  of.write("\nEND MODULE %s\n" % module_name)
