import tensorflow as tf
import numpy as np
import sys
import math


from binarization_utils import BinaryDense, AcidBaseDense, ph_activation

from keras.layers import Dense, Convolution2D, Activation, Flatten, MaxPooling2D,Input,Dropout,GlobalAveragePooling2D
from keras.layers.normalization import BatchNormalization
from binarization_utils import binarize_weight
from pathlib import Path



datatype = sys.argv[1]
assert datatype in ['bin', '3bit'] 
imgsize = int(sys.argv[2])
assert(imgsize in [8, 12, 16, 28])

train = len(sys.argv) >= 4





batch_norm_eps=1e-4
batch_norm_alpha=0.1#(this is same as momentum)
IMG_ROWS, IMG_COLS = imgsize, imgsize
input_size = IMG_ROWS * IMG_COLS
num_classes = 10

if train:
	num_classes = int(sys.argv[3] or 2)

Path('datasets/dataset_' + datatype + '/img' + str(imgsize) + '/').mkdir(parents=True, exist_ok=True)


print(datatype, imgsize, num_classes)

batch_size = 100


model = tf.keras.models.Sequential()
model.add(tf.keras.layers.InputLayer((input_size,), batch_size=batch_size))
model.add(BinaryDense(n_in=input_size, n_out=num_classes, gamma=1.0, input_shape=[input_size],width=16.0,point=0.0, first_layer=True))
model.add(Activation(ph_activation))



# learning rate schedule
def step_decay(epoch):
	initial_lrate = 0.025
	drop = 0.5
	epochs_drop = 50.0
	lrate = initial_lrate * math.pow(drop, math.floor((1+epoch)/epochs_drop))
	return lrate
	

def tf_log10(self, x):
		numerator = tf.math.log(x)
		denominator = tf.math.log(tf.constant(10, dtype=numerator.dtype))
		return numerator / denominator


def tohex(val, nbits):
  val = int(val * 255)
  return hex((val + (1 << nbits)) % (1 << nbits))[2:]


def filter_classes(inputs, labels, num_classes):
    inputs = np.array([value for index, value in enumerate(inputs) if labels[index] < num_classes], dtype=np.uint8)
    labels = np.array([label for label in labels if label < num_classes], dtype=np.uint8)
    return inputs, labels




mnist = tf.keras.datasets.mnist
(x_train, y_train), (x_test, y_test) = mnist.load_data()

if train:
	(x_train, y_train) = filter_classes(x_train, y_train, num_classes)
	(x_test, y_test) = filter_classes(x_test, y_test, num_classes)

new_batch_sz = x_train.shape[0] - (x_train.shape[0] % batch_size)
x_train, y_train = x_train[: new_batch_sz, :, :], y_train[: new_batch_sz] 
x_test, y_test = x_test[: new_batch_sz, :, :], y_test[: new_batch_sz] 

x_train = x_train.astype(np.float32)
x_test = x_test.astype(np.float32)


x_train = x_train.reshape(x_train.shape[0], 28, 28, 1)
x_test = x_test.reshape(x_test.shape[0], 28, 28, 1)

x_train = tf.image.resize_with_pad(x_train, IMG_ROWS, IMG_COLS).numpy()
x_test = tf.image.resize_with_pad(x_test, IMG_ROWS, IMG_COLS).numpy()

x_train = x_train.reshape(x_train.shape[0], IMG_ROWS, IMG_COLS)
x_test = x_test.reshape(x_test.shape[0], IMG_ROWS, IMG_COLS)


x_train, x_test = x_train / 255.0, x_test / 255.0 # normalizing data


x_train = x_train.reshape(-1, IMG_ROWS * IMG_COLS)
x_test = x_test.reshape(-1, IMG_ROWS * IMG_COLS)


y_train = tf.keras.utils.to_categorical(y_train, num_classes)
y_test = tf.keras.utils.to_categorical(y_test, num_classes)
discrete_layer = None

if datatype == '3bit':
	discrete_layer = tf.keras.layers.experimental.preprocessing.Discretization(bin_boundaries=[0.0, 1.0/7.0, 2.0/7.0, 3.0/7.0, 3.0/7.0 + 1e-10, 4.0/7.0, 5.0/7.0, 6.0/7.0, 7.0/7.0])
elif datatype == "bin":
	discrete_layer = tf.keras.layers.experimental.preprocessing.Discretization(bin_boundaries=[0.0])


if datatype == 'bin':
	x_train = 2 * discrete_layer(x_train) - 1
	x_test = 2 * discrete_layer(x_test) - 1
elif datatype == "3bit":
	x_train = discrete_layer(x_train) - 4
	x_test = discrete_layer(x_test) - 4


tf.print(tf.reshape(x_train[162], shape=(IMG_ROWS, IMG_COLS)), summarize=8)


num_map = {
	0: 'zero',
	1: 'one',
	2: 'two',
	3: 'three',
	4: 'four',
	5: 'five',
	6: 'six',
	7: 'seven',
	8: 'eight',
	9: 'nine'
}
counters = {
}
for digit, _ in num_map.items():
	counters[digit] = 1

lr = 0.0001
epochs = 50

# Choosing our optimizer
if train:
	optimizer = tf.keras.optimizers.Adam(learning_rate=lr)
	model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
	lrate = tf.keras.callbacks.ReduceLROnPlateau(
				monitor='val_accuracy', factor=0.5, patience=50, verbose=0, mode='auto',
				min_delta=0.0001, cooldown=0, min_lr=0
			)

	history_cb = model.fit(x_train, y_train, batch_size=batch_size, validation_data=(x_test, y_test), verbose=2,epochs=epochs,callbacks=[lrate])
	history = history_cb.history

	with open('datasets/dataset_' + datatype + '/kernel_acc_' + str(num_classes) + '_' + str(imgsize) + '.txt', mode='w') as acc_file:
		acc_file.write('loss:' + str(history['loss'][-1]) + ' accuracy:' + str(history['accuracy'][-1]) + ' val_loss:' + str(history['val_loss'][-1]) + ' val_accuracy:' + str(history['val_accuracy'][-1]))
		acc_file.write('\n')

	indeterm_count_p = 0
	indeterm_count_n = 0
	with open('datasets/dataset_' + datatype + '/kernel_mem_' + str(num_classes) + '_' + str(imgsize) + '.txt', mode='w') as weights_file:
		bin_weights = binarize_weight(model.layers[0].w)
		for c in range(num_classes):
			for ws in bin_weights:
				weights_file.write(str(int(ws[c])))
				weights_file.write('\n')
		with open('datasets/dataset_' + datatype + '/neurons_' + str(num_classes) + '_' + str(imgsize) + '.txt', mode='w') as neurons_file:
			with open('datasets/dataset_' + datatype + '/labels_' + str(num_classes) + '_' + str(imgsize) + '.txt', mode='w') as labels_file:
				for i, img in enumerate(x_test):
					result = tf.matmul([img], tf.cast(bin_weights, dtype=tf.int32))[0]
					lbl = int(tf.argmax(y_test[i]))
					tf_lbl = int(tf.argmax(result))
					filename = 'img_' + num_map[lbl] + '_' + str(counters[lbl]) + '.txt'
					labels_file.write(filename)
					labels_file.write(',')
					labels_file.write(str(tf_lbl))
					labels_file.write('\n')
					counters[lbl] = counters[lbl] + 1
					count_neg = 0
					count_pos = 0
					for r in result:
						if r <= 0:
							count_neg = count_neg + 1
						else:
							count_pos = count_pos + 1
					if count_neg != (num_classes - 1):
						indeterm_count_n = indeterm_count_n + 1

					neurons_file.write(filename)
					neurons_file.write(',')
					for j, neu in enumerate(result):
						neurons_file.write(str(float(neu)))
						if j != len(result) - 1:
							neurons_file.write(",")
					neurons_file.write('\n')
	tf.print(indeterm_count_p, indeterm_count_n, '/', len(x_test))


for i, img in enumerate(x_test):
	lbl = int(tf.argmax(y_test[i]))
	filename = 'datasets/dataset_' + datatype + '/img' + str(imgsize) + '/img_' + num_map[lbl] + '_' + str(counters[lbl]) + '.txt'
	if (i % 100) == 0:
		print(filename, i + 1, '/', len(x_test))
	counters[lbl] = counters[lbl] + 1
	with open(filename, mode='w') as imgfile:
		for px in img:
			imgfile.write(str(int(px)))
			imgfile.write('\n')
				