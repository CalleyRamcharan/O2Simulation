import numpy as np;
import scipy as sp;
import tensorflow as tf;
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Conv2D, MaxPooling2D
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import tensorflow.keras.callbacks as callbacks
from tensorflow.keras.callbacks import TensorBoard


def genArrays(filePath,index):
    images = []
    labels = []
    file = open(filePath,'r')

    header = file.readline()

    adcmax = []
    adchi= []
    adclo = []
    adchineighbour = []

    count2 = 0
    for line in file:
        
        if line == '\n':
            continue
        else:
            line = line.split(",")
            adcmax.append(int(line[0])/1023)
            adchi.append(int(line[1])/1023)
            adclo.append(int(line[2])/1023)
            adchineighbour.append(int(line[3])/1023)
            
            if(len(adcmax)==30):
                labels.append(index)
                images.append([])
                #newarr = [adcmax,adchi,adclo,adchineighbour]
                #images[count2].append(newarr)
                #images.append(newarr)
                images[count2].append(adcmax)
                images[count2].append(adchi)
                images[count2].append(adclo)
                images[count2].append(adchineighbour)
                adcmax = []
                adclo = []
                adchi = []
                adchineighbour = []
                count2 +=1
    return images, labels

#imgs = genArrays("./pion_csvs/testfile.csv")

#imgs = np.asarray(imgs)

#print((imgs.shape[0],imgs.shape[1],imgs.shape[2],1))

'''
for i in range(0,len(imgs)):
    for j in range(0,len(imgs[i])):
        print(imgs[i][j])
    print()'''

categories = [0,1] #zero is pion, one is electron

#generate features and labels
pion_imgs, pion_labels = genArrays("./pion_csvs/0049.csv",0)
p_0050, p_0050_l = genArrays("./pion_csvs/0050.csv",0)
p_0051, p_0051_l = genArrays("./pion_csvs/0051.csv",0)

pion_imgs = pion_imgs + p_0050 + p_0051
pion_labels = pion_labels + p_0050_l + p_0051_l

electron_imgs, electron_labels = genArrays("./electron_csvs/0049.csv",1)
e_0050, e_0050_l = genArrays("./electron_csvs/0050.csv",0)
e_0051, e_0051_l = genArrays("./electron_csvs/0051.csv",0)

electron_imgs = electron_imgs + e_0050 + e_0051
electron_labels = electron_labels + e_0050_l + e_0051_l

#ensuring 50/50 data fed into NN - more electron images than pion images so slice
electron_imgs = electron_imgs[:len(pion_imgs)]
electron_labels = electron_labels[:len(pion_imgs)]

#join the two arrays together
imgs = pion_imgs + electron_imgs
labels = pion_labels+electron_labels

#need to reshuffle data - generate random indices and reshuffle accordingly
imgs, labels = shuffle(imgs, labels, random_state=0)
imgs = np.array(imgs).reshape(-1,4,30,1)
labels = np.array(labels)

print("data loaded")
'''
pickle_out = open("imgs.pickle","wb")
pickle.dump(imgs,pickle_out)
pickle_out.close()

pickle_out = open("labels.pickle","wb")
pickle.dump(labels,pickle_out)
pickle_out.close()'''

model = Sequential()
model.add(Conv2D(64,(3,3), input_shape = imgs.shape[1:]))
model.add(Activation("relu"))
model.add(MaxPooling2D(pool_size=(2,2)))

'''
model.add(Conv2D(64,(3,3)))
model.add(Activation("relu"))
model.add(MaxPooling2D(pool_size=(2,2)))
'''
model.add(Flatten())

model.add(Dense(64))
model.add(Activation('relu'))

model.add(Dense(1))
model.add(Activation("sigmoid"))

model.compile(loss="binary_crossentropy",optimizer="adam",metrics=["accuracy"])

NAME = "e-vs-p-CNN2"
tensorboard = TensorBoard(log_dir="logs/{}".format(NAME))
earlystopping = callbacks.EarlyStopping(monitor ="val_loss", mode ="min", patience = 5, restore_best_weights = True)

model.fit(imgs,labels,batch_size=32,epochs=25,validation_split=0.3,callbacks=[tensorboard,earlystopping])

