#!/bin/bash

# Copyright (C) 2016 Stephan Kuschel <Stephan.Kuschel@gmail.com>
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

# THIS FILE MUST RUN WITHOUT ERROR ON EVERY COMMIT!

set -e

# chdir to epoch top level path
curdir=$(dirname "$(readlink -f "$0")")
cd $curdir/..

# Build SDF/C and install python sdf reader
(cd SDF/C; make)
SDF/utilities/build

# show system info
scripts/system_info.sh

# run the actual tests
scripts/run-tests.py 1d $@
scripts/run-tests.py 2d $@
scripts/run-tests.py 3d $@
