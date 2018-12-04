import os

def getfilenameindex(filename):
  filenum = 0
  while (os.path.exists(filename+"{:04}".format(filenum)+".png")):
    filenum+=1
  filename = "{}{:04}.png".format(filename,filenum)
  return filename

filename = "output/stream_"
filename = getfilenameindex(filename)
print(filename)
