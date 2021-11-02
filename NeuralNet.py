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

sim_images, sim_labels = genArrays("./muon_gun_10GeV/muonimages.csv",0)

cosmics_images, cosmics_labels = genArrays("./Cosmics/cosmicsimages.csv",1)

sim_images = sim_images[:35168]
sim_labels = sim_labels[:35168]

test_sim_images = sim_images[:50]
test_sim_labels = sim_labels[:50]

test_cosmics_images = cosmics_images[:50]
test_cosmics_labels = cosmics_labels[:50]

train_sim_images = sim_images[50:]
train_sim_labels = sim_labels[50:]

train_cosmics_images = cosmics_images[50:]
train_cosmics_labels = cosmics_labels[50:]

imgs = train_sim_images + train_cosmics_images
labels = train_sim_labels + train_cosmics_labels

imgs_test = test_cosmics_images + test_sim_images
labels_test = test_cosmics_labels + test_sim_labels

#need to reshuffle data - generate random indices and reshuffle accordingly
imgs, labels = shuffle(imgs, labels, random_state=0)
imgs = np.array(imgs).reshape(-1,4,30,1)
labels = np.array(labels)

imgs_test, labels_test = shuffle(imgs_test, labels_test, random_state=0)
imgs_test = np.array(imgs_test).reshape(-1,4,30,1)
labels_test = np.array(labels_test)

print("data loaded")

model = Sequential()
model.add(Conv2D(64,(3,3), input_shape = imgs.shape[1:]))
model.add(Activation("relu"))
model.add(MaxPooling2D(pool_size=(2,2)))

model.add(Flatten())

model.add(Dense(64))
model.add(Activation('relu'))

model.add(Dense(1))
model.add(Activation("sigmoid"))

model.compile(loss="binary_crossentropy",optimizer="adam",metrics=["accuracy"])

NAME = "cosmics-vs-sim"
tensorboard = TensorBoard(log_dir="logs/{}".format(NAME))
earlystopping = callbacks.EarlyStopping(monitor ="val_loss", mode ="min", patience = 5, restore_best_weights = True)

model.fit(imgs,labels,batch_size=32,epochs=3,validation_split=0.3,callbacks=[tensorboard,earlystopping])

prediction = model.predict(imgs_test)

for i in range(0,len(prediction)):
    print("pred = ",int(np.round(prediction[i][0],1)),"\t true = ",labels_test[i])