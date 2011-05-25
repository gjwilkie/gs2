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

def initialiseGlobalDimensions(outFile, inFile):
  for dim in inFile.dimensions:
    size = len(inFile.dimensions[dim])
    if(dim == splitDim):
      size = sum(splitSizes[dim])
    outFile.createDimension(dim , size)

def initialiseLocalDimensions(outFile, inFile, proc):
  for dim in inFile.dimensions:
    size = len(inFile.dimensions[dim])
    if(dim == splitDim):
      size = splitSizes[dim][proc]
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

def copyExtractSection(inFile, outFile):
  # Open the input file
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
        a = inVar[start:end,:,:]
        outVar[:] = a


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

def sizeSplit(inFile, nproc):

  for splitDim in splitDims:
    if(inFile.dimensions.has_key(splitDim)):
      worldLength = len(inFile.dimensions[splitDim])
      blocksize = (worldLength / nproc) + 1
      remainder = worldLength - (blocksize * (nproc-1))
      print remainder, blocksize
      for proc in xrange(0,nproc-1):
        splitSizes[splitDim].append(blocksize)
      splitSizes[splitDim].append(remainder)



   



def copyDecomposed(outFile, filePrefix, nprocs):
  # Write each of the identical files
  for proc in xrange(0,nprocs-1):
    copyDecomposedSection(outFile, filePrefix, proc)

if __name__ == "__main__":

  usage = "usage: %prog <nprocs> <fileprefix>"
  parser = optparse.OptionParser(usage)

  (options, args) = parser.parse_args()
  if(len(args) != 2):
    parser.error("Two arguments, the filename prefix and the number of processors are required")

  nprocs = int(args[0])
  filePrefix = args[1]
   
  inFile = netCDF4.Dataset(filePrefix + '.joined', 'r', format='NETCDF4')

  sizeSplit(inFile, nprocs)

  for proc in xrange(0,nprocs):
    filename = getFilename(filePrefix, proc)
    outFile = netCDF4.Dataset(filename, 'w', format='NETCDF3_CLASSIC')
    initialiseLocalDimensions(outFile, inFile, proc)
    initialiseVariables(outFile, inFile)
    copyRepeated(outFile, inFile)
    copyExtractSection(inFile, outFile)
    outFile.close()

#    copyDecomposed(outFile, filePrefix)

  inFile.close()
#  outFile.close()


