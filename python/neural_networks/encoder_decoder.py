import tensorflow as tf
from tensorflow.python.ops import control_flow_util
control_flow_util.ENABLE_CONTROL_FLOW_V2 = True
from tensorflow.keras.layers import MaxPooling2D, Concatenate, Conv2DTranspose, \
    Conv2D, Softmax, BatchNormalization, Lambda
from tensorflow.keras.models import Model
from .util.configurations import BatchIteratorConfiguration
from .data_io.label_rules import LabelRules
from .custom_layers import Conv2DNormRelu, Conv2DTransposeNormRelu
from .custom_losses import sample_weighted_loss_function


class EncoderDecoder(Model):
    """
    base line encoder decoder model. other models shouldn't do worse than this
    in semantic segmentation
    """
    def __init__(self, label_rules, n_pools=3,
                 n_start_filters=16, kernel_sizes=3, pool_sizes=2):
        super(EncoderDecoder, self).__init__()
        assert isinstance(label_rules, LabelRules)
        self.label_rules = label_rules
        self.n_labels = len(self.label_rules.input_colors)
        assert n_pools >= 1
        self.n_pools = n_pools
        assert n_start_filters >= 1
        self.n_start_filters = n_start_filters
        if isinstance(kernel_sizes, list):
            assert len(kernel_sizes) == n_pools
            self.kernal_sizes = kernel_sizes
        elif isinstance(kernel_sizes, int):
            self.kernal_sizes = [kernel_sizes] * n_pools
        else:
            raise ValueError('kernel_sizes should be a list or an integer')
        if isinstance(pool_sizes, list):
            assert len(pool_sizes) == n_pools
            self.pool_sizes = pool_sizes
        elif isinstance(pool_sizes, int):
            self.pool_sizes = [pool_sizes] * n_pools
        else:
            raise ValueError('kernel_sizes should be a list or an integer')

        self._build_model()

    def __str__(self):
        return self.__class__.__name__

    def _build_model(self):
        assert self.n_labels is not None
        n_filters = {'contract0': self.n_start_filters}
        for i in range(self.n_pools):
            if i > 0:
                n_filters['contract{}'.format(i)] = \
                    n_filters['contract{}'.format(i - 1)] * 2
            setattr(self, 'conv{0}_{1}x{1}s1'.format(i, self.kernal_sizes[i]),
                    Conv2DNormRelu(n_filters['contract{}'.format(i)], self.kernal_sizes[i], strides=1,
                                   negative_slope=0.1, name='conv2d{}'.format(i)))
            setattr(self, 'pool{0}_{1}x{1}s2'.format(i, self.pool_sizes[i]),
                    MaxPooling2D(pool_size=self.pool_sizes[i], strides=2,
                                 padding='same', data_format='channels_last', name='pool{}'.format(i)))
        for i in range(self.n_pools - 1, -1, -1):
            if i == self.n_pools - 1:
                n_filters['expand{}'.format(i)] = \
                    n_filters['contract{}'.format(i)]
            else:
                n_filters['expand{}'.format(i)] = \
                    n_filters['expand{}'.format(i + 1)] // 2
            setattr(self, 'convt{0}_{1}x{1}s2'.format(i, self.kernal_sizes[i]),
                    Conv2DTransposeNormRelu(n_filters['expand{}'.format(i)], self.kernal_sizes[i], strides=2,
                                            negative_slope=0.1, name='conv2dt{}'.format(i)))
        self.fc = Conv2D(self.n_labels, 1, strides=1, padding='same', data_format='channels_last', use_bias=False,
                         kernel_initializer=tf.initializers.GlorotNormal(), name='fc')
        self.fc_bn = BatchNormalization(axis=-1)
        self.softmax_activation = Softmax(axis=-1)
        self.classifier = Lambda(lambda x: tf.math.argmax(x, axis=-1, output_type=tf.int32))

    def call(self, inputs, training=None):
        for i in range(self.n_pools):
            if i == 0:
                result = getattr(self, 'conv{0}_{1}x{1}s1'.format(i, self.kernal_sizes[i]))(inputs, training=training)
            else:
                result = getattr(self, 'conv{0}_{1}x{1}s1'.format(i, self.kernal_sizes[i]))(result, training=training)
            result = getattr(self, 'pool{0}_{1}x{1}s2'.format(i, self.pool_sizes[i]))(result)
        for i in range(self.n_pools - 1, -1, -1):
            result = getattr(self, 'convt{0}_{1}x{1}s2'.format(i, self.kernal_sizes[i]))(result, training=training)
        result = self.fc(result)
        #result = self.fc_bn(result, training=training)
        result = self.softmax_activation(result)
        return result

    def compile_model(self, optimizer):
        self.compile(optimizer=optimizer,
                     loss=sample_weighted_loss_function(tf.losses.SparseCategoricalCrossentropy))

    def predict_from_inputs(self, inputs):
        activation = self.predict_on_batch(inputs)
        return self.predict_from_activation(activation)

    def predict_from_activation(self, activation):
        return self.classifier(activation)
