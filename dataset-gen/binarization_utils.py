# Original code at https://github.com/awai54st/BinaryNet-with-Keras-and-TF
import numpy as np
import pickle
import functools

import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Convolution2D, Activation, Flatten, MaxPooling2D,Input,Dropout,GlobalAveragePooling2D,Cropping2D,BatchNormalization
from tensorflow.keras.layers import Layer
#from multi_gpu import make_parallel


def FXP_quantize(
        value, point, width):
    point = tf.cast(point, 'float32')
    width = tf.cast(width, 'float32')
    maximum_value = tf.cast(2 ** (width - 1), 'float32')

    # x << (width - point)
    shift = 2.0 ** (tf.round(width) - tf.round(point))
    value_shifted = value * shift
    # quantize
    value_shifted = tf.round(value_shifted)
    value_shifted = tf.clip_by_value(value_shifted, -maximum_value, maximum_value - 1)
    # revert bit-shift earlier
    return value_shifted / shift, point, width

def log2(x):
    numerator = tf.log(x)
    denominator = tf.log(tf.constant(2, dtype=numerator.dtype))
    return numerator / denominator

def LOG_quantize(
        value, width):

    sign = tf.sign(value)
    value = log2(tf.abs(value))
    # quantize
    value, __, __ = FXP_quantize(
        value, width, width)
    # represent
    zero = tf.constant(0, dtype=tf.float32)
    return sign * (2 ** value)

def binarize_weight(x):
    '''Element-wise rounding to the closest integer with full gradient propagation.
    A trick from [Sergey Ioffe](http://stackoverflow.com/a/36480182)
    '''
    clipped = tf.clip_by_value(x, -1, 1)
    rounded = tf.sign(clipped)
    return clipped + tf.stop_gradient(rounded - clipped)

def binarize_activation(x):
    '''Element-wise rounding to the closest integer with full gradient propagation.
    A trick from [Sergey Ioffe](http://stackoverflow.com/a/36480182)
    '''
    clipped = tf.clip_by_value(x, -1, 1)
    rounded = tf.sign(clipped)
    return clipped + tf.stop_gradient(rounded - clipped)

def tf_custom_gradient_method(f):
    @functools.wraps(f)
    def wrapped(self, *args, **kwargs):
        if not hasattr(self, '_tf_custom_gradient_wrappers'):
            self._tf_custom_gradient_wrappers = {}
        if f not in self._tf_custom_gradient_wrappers:
            self._tf_custom_gradient_wrappers[f] = tf.custom_gradient(lambda *a, **kw: f(self, *a, **kw))
        return self._tf_custom_gradient_wrappers[f](*args, **kwargs)
    return wrapped


class Quantise_grad_layer(Layer):
    def __init__(self,width,point):
        super(Quantise_grad_layer, self).__init__()
        self.width = width
        self.point = point

    def build(self, input_shape):
        super(Quantise_grad_layer, self).build(input_shape)  # Be sure to call this at the end
        beta = np.zeros([5])*1.0
        self.beta=tf.variable(beta)

    def call(self, x):
        return self.quantise_gradient_op(x)

    @tf_custom_gradient_method
    def quantise_gradient_op(self, x):
        result = x # do forward computation
        def custom_grad(dy):
            grad, self.point, self.width = FXP_quantize(dy, self.point, self.width)
            return grad
        return result, custom_grad


class ResidualSign(Layer):
    def __init__(self, gamma=1,**kwargs):
        self.width = 16.0
        self.point = -2.0
        self.gamma=gamma
        super(ResidualSign, self).__init__(**kwargs)
    def call(self, x, mask=None):
        out_bin = binarize_activation(x)*tf.abs(self.gamma)
        return out_bin

    def get_output_shape_for(self,input_shape):
        return input_shape
    def compute_output_shape(self,input_shape):
        return input_shape


def switch(condition, t, e):
	import tensorflow as tf
	return tf.where(condition, t, e)


class BinaryConv(Layer):
	def __init__(self,nfilters,ch_in,k,padding,gamma,width,point,strides=(1,1),first_layer=False,**kwargs):
		self.nfilters=nfilters
		self.ch_in=ch_in
		self.k=k
		self.padding=padding
		self.strides=strides
		self.width = width
		self.point = point
		self.gamma=gamma
		self.first_layer = first_layer
		super(binary_conv,self).__init__(**kwargs)
	def build(self, input_shape):
		stdv=1/np.sqrt(self.k*self.k*self.ch_in)
		w = np.random.normal(loc=0.0, scale=stdv,size=[self.k,self.k,self.ch_in,self.nfilters]).astype(np.float32)
		if keras.backend._backend=="mxnet":
			w=w.transpose(3,2,0,1)
		self.w=tf.variable(w)
		self.trainable_weights=[self.w]


	def call(self, x,mask=None):
		self.out=tf.nn.conv2d(x, kernel=binarize_weight(self.w), padding=self.padding,strides=self.strides )
		#self.out = self.xnor_wg_conv_op(x,binarize_weight(self.w))
		self.output_dim=self.out.get_shape()
		return self.out

	@tf_custom_gradient_method
	def xnor_wg_conv_op(self, x, w):
		result=tf.nn.conv2d(x, kernel=w, padding=self.padding,strides=self.strides )
		def custom_grad(dy):

			w_reversed = tf.reverse(w, [0,1])
			w_reversed = tf.transpose(w_reversed, [0,1,3,2])

			# Valid
			dy_pad_number = tf.constant([[0, 0,], [2, 2,],[2, 2,], [0, 0]])
			dy_padded = tf.pad(dy, dy_pad_number)
			# if po2 dx
			#dy_padded = dy_padded * (2**(16))
			dy_padded = LOG_quantize(dy_padded, 8.0)
			#dy_padded = dy_padded * (2**(-16))
			dx = tf.nn.conv2d(dy_padded, kernel=w_reversed, padding="valid",strides=self.strides )

			if self.padding == "same":
				# Undo padding
				dx = Cropping2D(cropping=((1, 1), (1, 1)))(dx)
				# Pad x
				x_pad_number = tf.constant([[0, 0,], [1, 1,],[1, 1,], [0, 0]])
				x_trans = tf.transpose(tf.pad(x, x_pad_number), [3,1,2,0])
			elif self.padding == "valid":
				x_trans = tf.transpose(x, [3,1,2,0])
			dy_trans = tf.transpose(dy, [1,2,0,3])

			# Shift
			#dy_trans = dy_trans * (2**(16))
			dy_trans = LOG_quantize(dy_trans, 8.0)
			#dy_trans = dy_trans * (2**(-16))

			# Ternary dy
			ones = tf.ones_like(dy_trans)
			zeros = tf.zeros_like(dy_trans)

			if self.first_layer:
				dw = tf.nn.conv2d(x_trans, kernel=dy_trans, padding="valid",strides=self.strides )
			else:
				x_trans = tf.sign(x_trans + 1e-16) # Adding bias 1e-8 to sign function
				x_patches = tf.extract_image_patches(x_trans,
					[1, self.output_dim[1], self.output_dim[2], 1],
					[1, self.strides[0], self.strides[1], 1], [1, 1, 1, 1],
					padding="VALID")
				# CUDA impl
				#dw = stochastic_matmul_grad_gpu_module.stochastic_matmul_grad_gpu(tf.reshape(x_patches, [tf.shape(w)[2]*tf.shape(w)[0]*tf.shape(w)[1], tf.shape(dy)[0]*tf.shape(dy)[1]*tf.shape(dy)[2]]), tf.reshape(dy_trans, [-1, self.nfilters]))
				#dw = tf.reshape(dw, [tf.shape(w)[2], tf.shape(w)[0], tf.shape(w)[1], tf.shape(w)[3]])

				# Keras conv impl
				dw = tf.nn.conv2d(x_trans, kernel=dy_trans, padding="valid",strides=self.strides )
			dw = tf.transpose(dw, [1,2,0,3])

			N = self.k*self.k*self.ch_in*self.nfilters
			dw = 1./np.sqrt(N) * tf.sign(dw + 1e-16) # Adding bias 1e-8 to sign function
		
			return (dx, dw)
	
		return result, custom_grad


	def  get_output_shape_for(self,input_shape):
		return (input_shape[0], self.output_dim[1],self.output_dim[2],self.output_dim[3])
	def compute_output_shape(self,input_shape):
		return (input_shape[0], self.output_dim[1],self.output_dim[2],self.output_dim[3])

class BinaryDense(Layer):
	def __init__(self,n_in,n_out,gamma,width,point,first_layer=False,**kwargs):
		self.n_in = n_in
		self.n_out = n_out
		self.width = width
		self.point = point
		self.gamma = gamma
		self.first_layer = first_layer
		super(BinaryDense,self).__init__(**kwargs)

	def build(self, input_shape):
		# stdv=1 / np.sqrt(self.n_in)
		# w = np.random.normal(loc=0.0, scale=stdv,size=[self.n_in,self.n_out]).astype(np.float32)
		# self.w = tf.variable(w)
		self.w = self.add_weight(
            shape=(self.n_in, self.n_out),
            initializer="random_normal",
            trainable=True
        )

	def call(self, x, mask=None):
		# print("===")
		# print(self.n_in, self.n_out, self.w.shape, x.shape)

		#self.out = self.xnor_wg_dense_op(x,binarize_weight(self.w))
		self.out = tf.matmul(x, binarize_weight(self.w))
		# print(self.out.shape)
		# print("xxx")
		return self.out

	@tf_custom_gradient_method
	def xnor_wg_dense_op(self, x, w):
		result = tf.matmul(x, w)
		def custom_grad(dy):

			# Shift
			#dy = dy * (2**(16))
			dy = LOG_quantize(dy, 8.0)
			#dy = dy * (2**(-16))

			dx = tf.matmul(dy, tf.transpose(w))

			# Stochastic GPU impl
			if self.first_layer == True:
				#dw = stochastic_matmul_grad_gpu_module.stochastic_matmul_grad_gpu(tf.transpose(x), dy)
				dw = tf.matmul(tf.transpose(x), dy)
			else:
				#dw = stochastic_matmul_grad_gpu_module.stochastic_matmul_grad_gpu(tf.transpose(tf.sign(x + 1e-16)), dy)
				dw = tf.matmul(tf.transpose(tf.sign(x + 1e-16)), dy) # Adding bias 1e-8 to sign function
			
			N = self.n_in*self.n_out
			dw = 1./tf.sqrt(N) * tf.sign(dw + 1e-16) # Adding bias 1e-8 to sign function

			return (dx, dw)
	
		return result, custom_grad

	def  get_output_shape_for(self,input_shape):
		return (input_shape[0], self.n_out)
	def compute_output_shape(self,input_shape):
		return (input_shape[0], self.n_out)


@tf.function
def binary_activation(x, threshold=0):
	# result = tf.nn.sigmoid(x)
	result = tf.nn.softmax(x)
	# result = tf.nn.tanh(x)
	# result = tf.nn.softsign(x)
	# tf.print(x[0])
	# tf.print("xxxxxxxxxxxxx")
	# tf.print(result[0])
	# tf.print("=============")
	# if(x[0][0] > 0 and result[0][0] > 0 and x[0][1] > 0 and result[0][1] > 0):
	# 	tf.print("xxxxxxxxxxxxx")
	# 	tf.print(x[0])
	# 	tf.print(result[0])
	# 	tf.print("=============")
	# if(x[0][0] < 0 and result[0][0] < 0 and x[0][1] < 0 and result[0][1] < 0):
	# 	tf.print("xxxxxxxxxxxxx")
	# 	tf.print(x[0])
	# 	tf.print(result[0])
	# 	tf.print("=============")

	return result

def ph_activation(x):
	return binary_activation(x, 7)

class AcidBaseDense(Layer):
	def __init__(self,n_in,n_out,gamma,width,point,first_layer=False, transfer_volume=0.1, acid_ph=1.0, base_ph=13.0, water_ph=7.0, water_epsilon=1e-8,**kwargs):
		self.n_in = n_in
		self.n_out = n_out
		self.width = width
		self.point = point
		self.gamma = gamma
		self.first_layer = first_layer
		self.transfer_volume = transfer_volume
		self.acid_ph = acid_ph
		self.base_ph = base_ph
		self.water_ph = water_ph
		self.water_epsilon = water_epsilon
		super(AcidBaseDense, self).__init__(**kwargs)

	def build(self, input_shape):
		# stdv=1 / np.sqrt(self.n_in)
		# w = np.random.normal(loc=0.0, scale=stdv,size=[self.n_in,self.n_out]).astype(np.float32)
		# self.w = tf.variable(w)
		self.w = self.add_weight(
            shape=(self.n_in, self.n_out),
            initializer="random_normal",
            trainable=True,
			name="ph_weights"
        )
		
	def log10(self, x):
		numerator = tf.math.log(x)
		denominator = tf.math.log(tf.constant(10, dtype=numerator.dtype))
		return numerator / denominator


	def matmul_call(self, x, mask=None):
		self.out = tf.matmul(x, binarize_weight(self.w))
		return self.out


	def bad_call(self, x, mask=None):
		binary_weights = binarize_weight(self.w)


		batch_sz = x.shape[0]
		input_sz = x.shape[1]
		num_weights = binary_weights.shape[-1]
		bool_weights = tf.reshape(tf.less(binary_weights, 0), shape=(num_weights, input_sz))
		bool_weights = tf.broadcast_to(bool_weights, (batch_sz, num_weights, input_sz))
		
		x = tf.reshape(tf.tile(x, [1, num_weights]), [batch_sz, num_weights, input_sz])
		
		x = tf.where(bool_weights, x, 14.0 - x)

		is_acids = tf.less(x, self.water_ph)
		is_bases = tf.greater_equal(x, self.water_ph)

		acid_moles = tf.where(is_acids, (10 ** -x) * self.transfer_volume, tf.zeros(shape=is_acids.shape))
		base_moles = tf.where(is_bases, (1e-14 / (10 ** -x)) * self.transfer_volume, tf.zeros(shape=is_bases.shape))
		
		acid_moles_sum = tf.reduce_sum(acid_moles, axis=2)
		base_moles_sum = tf.reduce_sum(base_moles, axis=2)

		num_acid_units = tf.reduce_sum(tf.cast(is_acids, tf.int32), axis=2)
		num_base_units = tf.reduce_sum(tf.cast(is_bases, tf.int32), axis=2)
		acid_is_limiting = tf.less(acid_moles_sum, base_moles_sum)

		acid_limiting_h_conc = 1e-14/((base_moles_sum - acid_moles_sum) / (input_sz * self.transfer_volume))
		base_limiting_h_conc = (acid_moles_sum - base_moles_sum) / (input_sz * self.transfer_volume)

		h_conc = tf.where(acid_is_limiting, acid_limiting_h_conc, base_limiting_h_conc)
		resultant_ph = -self.log10(h_conc)

		self.out = resultant_ph
		return self.out

	def bad_call2(self, x, mask=None):
		binary_weights = binarize_weight(self.w)
		g = tf.identity(x)
		
		#x (100, 784)
		tf.print("ooooooo")
		tf.print(tf.round(tf.reshape(x[0], (28, 28))), summarize=28)



		batch_sz = x.shape[0]
		input_sz = x.shape[1]
		num_weights = binary_weights.shape[-1]
		print('binary_weights', binary_weights.shape)
		broadcasted_weights = tf.reshape(binary_weights, shape=(num_weights, input_sz))
		m = tf.reshape(self.w, shape=(num_weights, input_sz))
		print('broadcasted_weights.shape', broadcasted_weights.shape)
		broadcasted_weights = tf.broadcast_to(broadcasted_weights, (batch_sz, num_weights, input_sz))
		m = tf.broadcast_to(m, (batch_sz, num_weights, input_sz))
		print('broadcasted_weights.shape', broadcasted_weights.shape)
		print(x.shape)


		x = tf.reshape(tf.tile(x, [1, num_weights]), [batch_sz, num_weights, input_sz])
		tf.print(tf.reshape(broadcasted_weights[0][0], (28, 28)), summarize=28)
		tf.print(tf.reshape(m[0][0], (28, 28)), summarize=28)
		tf.print(tf.reshape((tf.math.multiply(x, broadcasted_weights) + 14.0)[0][0], (28, 28)), summarize=28)

		s = tf.reshape(tf.reduce_sum(tf.cast(broadcasted_weights, dtype=tf.float32), axis=0), (input_sz, num_weights)) 

		tf.print("yyyyyy")

		x = tf.math.mod(tf.math.multiply(x, broadcasted_weights) + 14.0, 14)

		tf.print(tf.round(tf.reshape(x[0][0], (28, 28))), summarize=28)
		tf.print("================")
		tf.print(tf.reshape(broadcasted_weights[0][0], (28, 28)), summarize=28)
		tf.print("===========")
		tf.print(tf.reshape(broadcasted_weights[0][5], (28, 28)), summarize=28)
		tf.print("===========")

		is_acids = tf.less(x, self.water_ph)
		is_bases = tf.greater_equal(x, self.water_ph)

		acid_moles = tf.where(is_acids, (10 ** -x) * self.transfer_volume, tf.zeros(shape=is_acids.shape))
		base_moles = tf.where(is_bases, (1e-14 / (10 ** -x)) * self.transfer_volume, tf.zeros(shape=is_bases.shape))
		
		acid_moles_sum = tf.reduce_sum(acid_moles, axis=2)
		base_moles_sum = tf.reduce_sum(base_moles, axis=2)

		acid_is_limiting = tf.less(acid_moles_sum, base_moles_sum)

		acid_limiting_h_conc = 1e-14/((base_moles_sum - acid_moles_sum) / (input_sz * self.transfer_volume))
		base_limiting_h_conc = (acid_moles_sum - base_moles_sum) / (input_sz * self.transfer_volume)

		h_conc = tf.where(acid_is_limiting, acid_limiting_h_conc, base_limiting_h_conc)
		resultant_ph = -self.log10(tf.abs(h_conc))
		

		tf.print(resultant_ph[0], summarize=10)
		tf.print(resultant_ph[1], summarize=10)
		tf.print((tf.matmul(g, binary_weights))[0], summarize=10)
		tf.print((tf.matmul(g, binary_weights))[1], summarize=10)
		tf.print("mmmmmmmmmmmmmmm")

		

		self.out = resultant_ph
		# print('resultant_ph.shape', resultant_ph.shape)
		return self.out

	def call(self, x, mask=None):
		binary_weights = binarize_weight(self.w)

		input_sz = x.shape[1]
		
		input_h_conc = 10**-x
		input_oh_conc = 1e-14/10**-x
		input_h_num_moles = input_h_conc * self.transfer_volume
		input_oh_num_moles = input_oh_conc * self.transfer_volume

		acid_moles_sum = tf.matmul(input_h_num_moles, binary_weights)
		base_moles_sum = tf.matmul(input_oh_num_moles, binary_weights)

		
# - Storage of Information using Small Organic Molecules 
# - Nature Communications 
# - ACS Central 


		remaining_moles = acid_moles_sum - base_moles_sum # -ve means acid is limiting agent, result is basic

		acid_is_limiting = tf.less(remaining_moles, 0)

		remaining_moles_conc = tf.abs(remaining_moles) / (input_sz * self.transfer_volume)

		tf.print("input", tf.round(tf.reshape(x[0], (28, 28))), summarize=28)
		tf.print(binary_weights, summarize=28)
		tf.print("remaining_moles_conc", remaining_moles_conc[0], summarize=10)

		resultant_h_conc = tf.where(acid_is_limiting, 1e-14/remaining_moles_conc, remaining_moles_conc)

		tf.print("resultant_h_conc", resultant_h_conc[0], summarize=10)
		
		resultant_ph = -self.log10(resultant_h_conc)

		tf.print("resultant_ph", resultant_ph[0], summarize=10)

		tf.print("==================")
		self.out = resultant_ph

		return self.out


	

	@tf_custom_gradient_method
	def xnor_wg_dense_op(self, x, w):
		result = tf.matmul(x, w)
		def custom_grad(dy):

			# Shift
			#dy = dy * (2**(16))
			dy = LOG_quantize(dy, 8.0)
			#dy = dy * (2**(-16))

			dx = tf.matmul(dy, tf.transpose(w))

			# Stochastic GPU impl
			if self.first_layer == True:
				#dw = stochastic_matmul_grad_gpu_module.stochastic_matmul_grad_gpu(tf.transpose(x), dy)
				dw = tf.matmul(tf.transpose(x), dy)
			else:
				#dw = stochastic_matmul_grad_gpu_module.stochastic_matmul_grad_gpu(tf.transpose(tf.sign(x + 1e-16)), dy)
				dw = tf.matmul(tf.transpose(tf.sign(x + 1e-16)), dy) # Adding bias 1e-8 to sign function
			
			N = self.n_in*self.n_out
			dw = 1./tf.sqrt(N) * tf.sign(dw + 1e-16) # Adding bias 1e-8 to sign function

			return (dx, dw)
	
		return result, custom_grad

	def  get_output_shape_for(self,input_shape):
		return (input_shape[0], self.n_out)
	def compute_output_shape(self,input_shape):
		return (input_shape[0], self.n_out)

class my_flat(Layer):
	def __init__(self,**kwargs):
		super(my_flat,self).__init__(**kwargs)
	def build(self, input_shape):
		return

	def call(self, x, mask=None):
		self.out=tf.reshape(x,[-1,np.prod(x.get_shape().as_list()[1:])])
		return self.out
	def  compute_output_shape(self,input_shape):
		shpe=(input_shape[0],int(np.prod(input_shape[1:])))
		return shpe