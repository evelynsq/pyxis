#!/bin/bash

# This script contains command line tests for pyxis

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

runcmd_pass()
{
    echo "[runcmd_pass]: $1"
    sh -c "$1" >/dev/null 2>&1 || die "Error running: $1"
}

runcmd_fail()
{
    echo "[runcmd_fail]: $1"
    sh -c "$1" >/dev/null 2>&1 && die "Command should have failed: $1"
}

if [ $# -eq 0 ]; then
    # use default example location
    EXDATADIR="example-files"
elif (( $# != 1 )) ; then
    echo "usage: cmdline_tests.sh {example_dir}" 2>&1
    echo "Expected 1 arguments but recieved $#" 2>&1
    exit 1
else
    EXDATADIR=$1
fi

TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

# Tests for pyxis

## Testing that basic usage passes with example files
runcmd_pass "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/test.pwms"
## Addition of background
runcmd_pass "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/test.pwms -b ${EXDATADIR}/background.bed"
## Addition of sequence logo generation
runcmd_pass "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/test.pwms -b ${EXDATADIR}/background.bed -s"

## Incorrect peaks file specified
runcmd_fail "pyxis ${EXDATADIR}/XYZA.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/test.pwms -b ${EXDATADIR}/background.bed"
## Incorrect reference genome file specified
runcmd_fail "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/XYZA.fa ${EXDATADIR}/test.pwms -b ${EXDATADIR}/background.bed"
## Incorrect PWMS file specified
runcmd_fail "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/XYZA.pwms -b ${EXDATADIR}/background.bed"
## Incorrect background peaks file specified
runcmd_fail "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa ${EXDATADIR}/test.pwms -b ${EXDATADIR}/XYZA.bed"
## One required positional argument not specified
runcmd_fail "pyxis ${EXDATADIR}/peaks.bed ${EXDATADIR}/ref.fa"
## Required positional arguments specified in incorrect order
runcmd_fail "pyxis ${EXDATADIR}/ref.fa ${EXDATADIR}/peaks.bed ${EXDATADIR}/test.pwms"
