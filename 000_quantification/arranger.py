import os, shutil

dir = '/hpcdata/Mimir/adrian/research/akthelia/data/'
files = os.listdir(dir)
print(files)

labels = []
for file in files:
    v = file.split('_')
    new = '_'.join(v[:4])
    labels.append(new)
uniquelabels = list(set(labels))
uniquelabels.sort()
print(len(uniquelabels), uniquelabels)

for uniquelabel in uniquelabels:
    new_path = dir + uniquelabel
    if os.path.exists(new_path) == False:
        os.mkdir(new_path)

    shutil.move('{}{}_1.fq.gz'.format(dir, uniquelabel), "{}/".format(new_path))
    shutil.move('{}{}_2.fq.gz'.format(dir, uniquelabel), "{}/".format(new_path))