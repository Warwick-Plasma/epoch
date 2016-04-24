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

# All arguments given to this script will be forwarded
# to 'epoch{1,2,3}d/run-tests.py'.

# Build SDF/C and install python sdf reader
(cd SDF/C; make)
SDF/utilities/build

# run the actual tests
epoch1d/run-tests.py $@
