#!/bin/sh
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_calib_262548_262799/* > inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_calib_262800_263230/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_calib_263231_263359/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_calib_263360_263379/* >> inlist.dat
ls -1 /rfs/sanders/crab_projects/crab_PbPb_2015_calib_263380_263614/* >> inlist.dat
root -l -b -q ClearBadFiles.C+

