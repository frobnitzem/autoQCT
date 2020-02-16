import re

###### Defining atoms ######
# All atoms will have a time step, i, part of a cluster energy, ener, and position <x, y, z>
class atom:
    def __init__(self,i,ener,type,x,y,z):
        self.i = i
        self.ener = ener
        self.type = type
        self.x = x
        self.y = y
        self.z = z

    def to_dict(self):
        return {
            'i': self.i,
            'ener': self.ener,
            'type': self.type,
            'x': self.x,
            'y': self.y,
            'z': self.z,
        }

    def __str__(self):
        return (str(self.i) + "\t" + str(self.ener) + "\t" +  str(self.type) + "\t" +  str(self.x) +  "\t" + str(self.y)
                + "\t" +  str(self.z))


# Reading data in
# This is a little hacky
def dataIn(cp2k, pattern):
    with open(cp2k) as file_in:
        maxFrame = 0
        dataFrameList = []
        frame = atom(0, 0, "NULL", 0, 0, 0)
        for line in file_in:
            for m in re.finditer(pattern, line):
                lines = line.split('\t')
                if re.search(r'i =', lines[0]):
                    temp = lines[0].split('=')
                    i = int(temp[1].split(',')[0])
                    ener = float(temp[3].split(',')[0])
                    clusterSize = 0
                else:
                    clusterSize = clusterSize + 1
                    temp = lines[0].split(" ")
                    temp = list(filter(None, temp))
                    frame.i = i
                    frame.ener = ener
                    frame.type = temp[0]
                    frame.x = float(temp[1])
                    frame.y = float(temp[2])
                    frame.z = float(temp[3])
                    appFrame = atom(frame.i, frame.ener, frame.type, frame.x, frame.y, frame.z)
                    dataFrameList.append(appFrame)
                    if frame.i > maxFrame:
                        maxFrame = frame.i

        return [dataFrameList, maxFrame, clusterSize]
