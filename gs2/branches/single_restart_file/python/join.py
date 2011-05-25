#!/usr/bin/env python

import netCDF4
import sys
import numpy
import os
import optparse

splitDims = ["glo"]
splitSizes = {};

def getFilename(filePrefix, proc):
  format = filePrefix + ".%d"
  return format % (proc)

for splitDim in splitDims:
  splitSizes[splitDim] = []

def initialiseDimensions(outFile, inFile):
  for dim in inFile.dimensions:
    size = len(inFile.dimensions[dim])
    if(dim == splitDim):
      size = sum(splitSizes[dim])
    outFile.createDimension(dim , size)

def initialiseVariables(outFile, inFile):
  for varName in inFile.variables:
    inVar = inFile.variables[varName]
    outFile.createVariable(varName, inVar.dtype, inVar.dimensions)


def copyRepeated(outFile, inFile):
  for varName in inFile.variables:
    inVar = inFile.variables[varName]
    outVar = outFile.variables[varName]
    if(not (splitDim in inVar.dimensions)):
      a = inVar[:]
      outVar[:] = a

def copyDecomposedSection(outFile, filePrefix, proc):
  # Open the input file
  filename = getFilename(filePrefix, proc)
  inFile = netCDF4.Dataset(filename, 'r', format='NETCDF4')

  for splitDim in splitDims:
    start = 0
    for x in xrange(0,proc):
      start += splitSizes[splitDim][x]
    end = start + splitSizes[splitDim][proc]
    print proc, start, end, end - start
    for varName in inFile.variables:
      outVar = outFile.variables[varName]
      inVar = inFile.variables[varName]
      if(splitDim in inVar.dimensions):
        a = inVar[:]
        outVar[start:end,:,:] = a
  inFile.close()

def scanFiles(filePrefix):

  proc = 0
  while True:
    filename = getFilename(filePrefix, proc)
    if(not os.path.isfile(filename)):
      break
    print "Scanning file " + filename
    inFile = netCDF4.Dataset(filename, 'r', format='NETCDF4')
    for splitDim in splitDims:
      if(inFile.dimensions.has_key(splitDim)):
        splitSizes[splitDim].append(len(inFile.dimensions[splitDim]))
    inFile.close()
    proc += 1
  return proc-1

def copyDecomposed(outFile, filePrefix):
  # Write each of the identical files
  for proc in xrange(0,nprocs+1):
    copyDecomposedSection(outFile, filePrefix, proc)



if __name__ == "__main__":

  usage = "usage: %prog <fileprefix>"
  parser = optparse.OptionParser(usage)

  (options, args) = parser.parse_args()
  if(len(args) != 1):
    parser.error("First argument, the filename prefix is required")

  filePrefix = args[0]
   

  nprocs = scanFiles(filePrefix)

  inFile = netCDF4.Dataset(getFilename(filePrefix, 0), 'r', format='NETCDF4')
  outFile = netCDF4.Dataset(filePrefix + '.joined', 'w', format='NETCDF4')

  initialiseDimensions(outFile, inFile)
  initialiseVariables(outFile, inFile)

  copyRepeated(outFile, inFile)
  copyDecomposed(outFile, filePrefix)

  inFile.close()
  outFile.close()


