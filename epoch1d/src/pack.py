#! /usr/bin/env python

# Copyright (C) 2009-2019 University of Warwick
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
import codecs
try:
    import argparse
    got_argparse = True
except:
    got_argparse = False


def str2bool(x):
    if x.lower() not in ['true', 'yes', '1', 'false', 'no', '0']:
        raise TypeError("Argument is not a Boolean string")
    return x.lower() in ['true', 'yes', '1']


def stripped(s):
    try:
        # python 2
        s = unicode(s, 'ascii', 'ignore')
    except(NameError):
        # python 3
        s = s.encode('ascii', 'ignore').decode()
    return s.strip()

if got_argparse:
    argp = argparse.ArgumentParser(
        description="Pack source code for writing to SDF output")
    argp.add_argument("prefix", type=str, help="Package name")
    argp.add_argument("pack_source_code", type=str2bool,
                      help="Pack source code")
    argp.add_argument("pack_git_diff", type=str2bool,
                      help="Pack git diff")
    argp.add_argument("pack_git_diff_from_origin", type=str2bool,
                      help="Pack git diff from origin")
    argp.add_argument("generate_checksum", type=str2bool,
                      help="Generate checksum")
    argp.add_argument("f77_output", type=str2bool,
                      help="Fortran 77 output")
    argp.add_argument("outfile", type=str, help="Output file")
    argp.add_argument("compiler_info", type=stripped, help="Compiler info")
    argp.add_argument("compiler_flags", type=stripped, help="Compiler flags")
    argp.add_argument("filelist", type=str, nargs='*', help="Source files")
    args = argp.parse_args()
else:
    args = type("", (), dict(dummy=1))()
    args.prefix = sys.argv[2]
    (args.pack_source_code,
     args.pack_git_diff,
     args.pack_git_diff_from_origin,
     args.generate_checksum,
     args.f77_output,) = map(str2bool, sys.argv[3:8])
    args.outfile = sys.argv[8]
    (args.compiler_info,
     args.compiler_flags,) = map(stripped, sys.argv[9:11])
    args.filelist = sys.argv[11:]

prefix = args.prefix
pack_source_code = args.pack_source_code
pack_git_diff = args.pack_git_diff
pack_git_diff_from_origin = args.pack_git_diff_from_origin
generate_checksum = args.generate_checksum

commitfile = os.path.join(os.environ['GIT_WORK_TREE'], 'src', 'COMMIT')

archive = "source_info_archive.tgz"
hexdump = "source_info_hexdump.txt"
gitdiff = "source_info_gitdiff.txt"
varname = "%s_bytes" % prefix
diffname = "%s_diff_bytes" % prefix
module_name = "%s_source_info" % prefix
outfile = args.outfile
incfile = os.path.splitext(outfile)[0] + '_include.inc'

f77_output = args.f77_output
nbytes = 8
nelements = 0
padding = 0
vname = varname

if f77_output:
    inc_handle = open(incfile, "w")
    linestart = 6*' '
    linecont = 5*' '+'&'
    suffix = ''
    ncolumns = 72
    ncontinuation = 19
else:
    linestart = ''
    linecont = ''
    suffix = '&'
    ncolumns = 132          # gfortran ignores the F90 standard of 139
    ncontinuation = 39


def byteswap4(s):
    s = binascii.unhexlify(s)
    a, = struct.unpack('>L', s)
    s = struct.pack('<L', a)
    return binascii.hexlify(s).decode('utf-8')


def byteswap8(s):
    s = binascii.unhexlify(s)
    a, = struct.unpack('>Q', s)
    s = struct.pack('<Q', a)
    return binascii.hexlify(s).decode('utf-8')


def byteswap(n, s):
    """Simplifies the logic of calling byteswap4/8"""

    if n == 4:
        return byteswap4(s)
    elif n == 8:
        return byteswap8(s)
    else:
        raise ValueError("Invalid value", nbytes)


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


def print_character(name, value):
    global vname, of
    # ilen = len(value)
    ilen = 256
    var = "%s_%s" % (vname, name)
    if f77_output:
        of = inc_handle
        wrapped("CHARACTER*%i %s" % (ilen, var))
        wrapped("COMMON/c_%s/%s" % (vname, var))
        of = out_handle
        wrapped("CHARACTER*%i %s" % (ilen, var))
        wrapped("COMMON/c_%s/%s" % (vname, var))
        wrapped("DATA %s/'%s'/" % (var, value))
    else:
        ilen = len(value)
        if ilen == 0:
            ilen = 1
        wrapped("CHARACTER(LEN=%i) :: %s = '%s'" % (ilen, var, value))


def print_integer(name, value):
    global vname, of
    var = "%s_%s" % (vname, name)
    if f77_output:
        of = inc_handle
        wrapped("INTEGER " + var)
        wrapped("COMMON/i_%s/%s" % (vname, var))
        of = out_handle
        wrapped("INTEGER " + var)
        wrapped("COMMON/i_%s/%s" % (vname, var))
        wrapped("DATA %s/%i/" % (var, value))
    else:
        wrapped("INTEGER, PARAMETER :: %s = %i" % (var, value))


def get_bytes_checksum(files):
    global checksum_type
    import codecs
    if not generate_checksum:
        checksum_type = ''
        return ''
    cksum = hashlib.new('sha256')
    for name in files:
        with codecs.open(name, encoding='utf-8') as f:
            while True:
                data = f.read(cksum.block_size)
                if not data:
                    break
                cksum.update(data.encode('utf-8'))
    checksum_type = 'sha256'
    return cksum.hexdigest()


def write_data_bytes(filename, varname):
    global mimetype, of
    global linestart, linecont, suffix, ncolumns, ncontinuation

    f = open(filename, 'rb')
    d = f.read()
    dhex = codecs.encode(d, 'hex_codec').decode('utf-8')
    f.close()
    os.remove(filename)

    nelements = (len(d)+nbytes-1) // nbytes
    padding = nelements * nbytes - len(d)
    dhex += '00' * padding

    print_character('mimetype', mimetype)
    print_integer('padding', padding)
    print_integer('len', nelements)
    print_integer_array(nelements)

    nwidth = len("z'',") + 2 * nbytes
    nper_line_body = (ncolumns - 1) // nwidth
    sdata = linestart + "DATA(%s(i),i=%i,%i)/" % (varname, nelements,
                                                  nelements)
    nper_line_first = (ncolumns - len(sdata) - 1) // nwidth
    nper_segment = nper_line_first + nper_line_body * ncontinuation

    i0 = 0
    segline = 0
    elements_written = 0
    while elements_written < nelements:
        ss = ""
        if segline == 0:
            i1 = min(i0 + nper_segment, nelements)
            ss += linestart + "DATA(%s(i),i=%i,%i)/" % (varname, i0+1, i1)
            i0 = i1
        else:
            ss += linecont

        shex = dhex[2*nbytes*elements_written:2*nbytes*(elements_written+1)]
        shex = byteswap(nbytes, shex)
        ss += "z'%s'" % shex
        elements_written = elements_written + 1
        if segline == 0:
            nper_line = nper_line_first - 1
        else:
            nper_line = nper_line_body - 1
        n = 0
        while n < nper_line and elements_written != nelements:
            shex = dhex[2*nbytes*elements_written:
                        2*nbytes*(elements_written+1)]
            shex = byteswap(nbytes, shex)
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
    global of, nbytes, vname
    if value == 0:
        value = 1
    if f77_output:
        of = inc_handle
        wrapped("INTEGER*%i %s(%s_len)" % (nbytes, vname, vname))
        of = out_handle
        wrapped("INTEGER*%i %s(%i)" % (nbytes, vname, value))
    else:
        wrapped("INTEGER(%i) :: %s(%i)" % (nbytes, vname, value))


try:
    cmd = sp.Popen("git describe --always --long --dirty", shell=True,
                   stderr=sp.PIPE, stdout=sp.PIPE)
    output = cmd.communicate()
    if cmd.returncode == 127:
        print('WARNING: Git command not found')
        git_version = ''
        pack_git_diff = False
        try:
            f = open(commitfile, "r")
            string = f.readline().rstrip('\n')
            f.close()
            git_version = string.split('=')[1]
        except:
            pass
    elif cmd.returncode != 0 and str(output[1]).find('ot a git repo') != -1:
        print('WARNING: Not a git repository')
        git_version = ''
        pack_git_diff = False
        try:
            f = open(commitfile, "r")
            string = f.readline().rstrip('\n')
            f.close()
            git_version = string.split('=')[1]
        except:
            pass
    elif cmd.returncode != 0:
        raise Exception('ERROR: unable to generate git diff')
    else:
        git_version = output[0].decode('utf-8').rstrip()
        pack_git_diff = True
except:
    raise Exception('ERROR: unable to generate git diff')

tsec = time.time()
compile_date = int(round(tsec))
compile_date_string = time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime(tsec))
compile_machine_info = ' '.join((platform.node(), platform.platform()))
compiler_info = args.compiler_info
compiler_flags = args.compiler_flags

# fnull=open(os.devnull,'w')
# filelist = sp.check_output("git ls-files --cached --no-empty-directory "
#                            + "--full-name", shell=True,
#                            stderr=fnull).rstrip()

filelist = args.filelist
if filelist == []:
    pack_source_code = False


# Now write the file

out_handle = open(outfile, "w")
of = out_handle

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
        sp.call(["git diff origin/main > %s" % gitdiff], shell=True)
    else:
        sp.call(["git diff > %s" % gitdiff], shell=True)
    if os.path.getsize(gitdiff) != 0:
        # Remove TABLES directory from diff
        with codecs.open(gitdiff, 'r', encoding='utf-8') as fd:
            lines = fd.readlines()
        with codecs.open(gitdiff, 'w', encoding='utf-8') as fd:
            ignore = False
            for l in lines:
                if l.startswith("diff"):
                    if l.find("/TABLES/") != -1:
                        ignore = True
                    else:
                        ignore = False
                if not ignore:
                    fd.write(l)

        checksum = get_bytes_checksum([gitdiff])
        zgitdiff = gitdiff + '.gz'
        with open(gitdiff, 'rb') as f_in:
            with gzip.open(zgitdiff, 'wb') as f_out:
                f_out.writelines(f_in)
        os.remove(gitdiff)
        os.rename(zgitdiff, gitdiff)
    mimetype = 'application/x-gzip'

    print_character('checksum_type', checksum_type)
    print_character('checksum', checksum)

    write_data_bytes(gitdiff, vname)

if f77_output:
    of.write(linestart + "END SUBROUTINE\n")
else:
    of.write("\nEND MODULE %s\n" % module_name)
