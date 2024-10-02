import tensorflow as tf
from tensorflow.python.ops import control_flow_util
control_flow_util.ENABLE_CONTROL_FLOW_V2 = True
from tensorflow.keras.layers import MaxPooling2D, Concatenate, Conv2DTranspose, \
    Conv2D, Softmax, BatchNormalization, Lambda, ReLU
from tensorflow.keras.models import Model
from .data_io.label_rules import LabelRules
from .custom_layers import Conv2DNormRelu, Conv2DTransposeNormRelu, Inception
from .custom_losses import sample_weighted_loss_function


class NisslNet(Model):
    def __init__(self, label_rules):
        super(NisslNet, self).__init__()
        assert isinstance(label_rules, LabelRules)
        self.label_rules = label_rules
        self.n_labels = len(self.label_rules.input_colors)

        self._build_model()

    def __str__(self):
        return self.__class__.__name__

    """
    def _predict_sigmoid(self):
        batch_activation = tf.nn.sigmoid(self.fc)
        predicted_labels = tf.round(batch_activation)
        predicted_labels = tf.cast(predicted_labels, dtype=tf.int32)
        return tf.reshape(predicted_labels, (self.batch_size, self.patch_height,
                                             self.patch_width))
                                             """

    def _build_model(self):
        assert self.n_labels is not None

        contract_filters0 = 32

        self.inception0 = Inception([contract_filters0 // 4 * 3, contract_filters0 // 4], [5, 7], reduced_dimensions=[24, 8],
                                    relu_negative_slope=0.1, name='inception0')

        self.pool_3x3s2 = MaxPooling2D(pool_size=(3, 3), strides=(2, 2), padding='same',
                                       data_format='channels_last', name='pool_3x3s2')
        self.concatenate = Concatenate(axis=-1)

        contract_filters1 = contract_filters0 * 2
        self.inception1 = Inception([contract_filters1 // 4 * 3, contract_filters1 // 4], [3, 5], reduced_dimensions=[48, 16],
                                    has_activation=False, name='inception1')
        self.res_conv2d1_3x3s2 = Conv2D(contract_filters1, 3, strides=2, padding='same', data_format='channels_last',
                                        use_bias=False, kernel_initializer=tf.initializers.GlorotNormal(), name='res_conv2d1')
        self.relu1 = ReLU(negative_slope=0.1)

        contract_filters2 = contract_filters1 * 2
        self.inception2 = Inception([contract_filters2 // 4 * 3, contract_filters2 // 4], [3, 5], reduced_dimensions=[96, 32],
                                    relu_negative_slope=0.1, name='inception2')

        contract_filters3 = contract_filters2 * 2
        self.inception3 = Inception([contract_filters3 // 4 * 3, contract_filters3 // 4], [3, 5], reduced_dimensions=[192, 64],
                                    has_activation=False, name='inception3')
        self.res_conv2d3_3x3s2 = Conv2D(contract_filters3, 3, strides=2, padding='same', data_format='channels_last',
                                        use_bias=False, kernel_initializer=tf.initializers.GlorotNormal(), name='res_conv2d3')
        self.relu3 = ReLU(negative_slope=0.1)

        contract_filters4 = contract_filters3 * 2
        self.inception4 = Inception([contract_filters4 // 4 * 3, contract_filters4 // 4], [3, 5], reduced_dimensions=[384, 128],
                                    relu_negative_slope=0.1, name='inception4')

        contract_filters5 = contract_filters4 * 3
        self.conv2d5_3x3s1 = Conv2DNormRelu(contract_filters5, 3, strides=1, negative_slope=0.1, name='conv2d5')

        expand_filters4 = contract_filters5 // 2
        self.conv2dt4_3x3s2 = Conv2DTransposeNormRelu(expand_filters4, 3, strides=2, negative_slope=0.1, name='conv2dt4')
        self.conv2d4_3x3s1 = Conv2DNormRelu(expand_filters4, 3, strides=1, negative_slope=0.1, name='conv2d4')

        expand_filters3 = expand_filters4 // 2
        self.conv2dt3_3x3s2 = Conv2DTransposeNormRelu(expand_filters3, 3, strides=2, negative_slope=0.1, name='conv2dt3')
        self.conv2d3_3x3s1 = Conv2DNormRelu(expand_filters3, 3, strides=1, negative_slope=0.1, name='conv2d3')

        expand_filters2 = expand_filters3 // 2
        self.conv2dt2_3x3s2 = Conv2DTransposeNormRelu(expand_filters2, 3, strides=2, negative_slope=0.1, name='conv2dt2')
        self.conv2d2_3x3s1 = Conv2DNormRelu(expand_filters2, 3, strides=1, negative_slope=0.1, name='conv2d2')

        expand_filters1 = expand_filters2 // 2
        self.conv2dt1_3x3s2 = Conv2DTransposeNormRelu(expand_filters1, 3, strides=2, negative_slope=0.1, name='conv2dt1')
        self.conv2d1_3x3s1 = Conv2DNormRelu(expand_filters1, 3, strides=1, negative_slope=0.1, name='conv2d1')

        expand_filters0 = expand_filters1 // 2
        self.conv2dt0_3x3s2 = Conv2DTransposeNormRelu(expand_filters0, 3, strides=2, negative_slope=0.1, name='conv2dt0')
        self.conv2d0_3x3s1 = Conv2DNormRelu(expand_filters0, 3, strides=1, negative_slope=0.1, name='conv2d0')

        self.fc = Conv2D(self.n_labels, 1, strides=1, padding='same', data_format='channels_last',
                         use_bias=False, kernel_initializer=tf.initializers.GlorotNormal(), name='fc')
        self.fc_bn = BatchNormalization(axis=-1)
        self.softmax_activation = Softmax(axis=-1)
        self.classifier = Lambda(lambda x: tf.math.argmax(x, axis=-1, output_type=tf.int32))

    def call(self, inputs, training=None):
        inception0 = self.inception0(inputs, training=training)

        pool0 = self.pool_3x3s2(inception0)
        inception1 = self.inception1(pool0, training=training)

        res1 = self.res_conv2d1_3x3s2(inputs) + inception1

        pool1 = self.pool_3x3s2(res1)
        inception2 = self.inception2(pool1, training=training)

        pool2 = self.pool_3x3s2(inception2)
        inception3 = self.inception3(pool2, training=training)

        res3 = self.res_conv2d3_3x3s2(pool1) + inception3

        pool3 = self.pool_3x3s2(res3)
        #inception4 = self.inception4(pool3, training=training)

        #pool4 = self.pool_3x3s2(inception4)
        #conv2dt4 = self.conv2dt4_3x3s2(pool4)
        #concat4 = self.concatenate([conv2dt4, inception4])
        #conv2d4 = self.conv2d4_3x3s1(concat4, training=training)

        #conv2dt3 = self.conv2dt3_3x3s2(conv2d4)
        conv2dt3 = self.conv2dt3_3x3s2(pool3)
        concat3 = self.concatenate([conv2dt3, inception3])
        conv2d3 = self.conv2d3_3x3s1(concat3, training=training)

        conv2dt2 = self.conv2dt2_3x3s2(conv2d3)
        concat2 = self.concatenate([conv2dt2, inception2])
        conv2d2 = self.conv2d2_3x3s1(concat2, training=training)

        conv2dt1 = self.conv2dt1_3x3s2(conv2d2)
        concat1 = self.concatenate([conv2dt1, inception1])
        conv2d1 = self.conv2d1_3x3s1(concat1, training=training)

        conv2dt0 = self.conv2dt0_3x3s2(conv2d1)
        concat0 = self.concatenate([conv2dt0, inception0])
        conv2d0 = self.conv2d0_3x3s1(concat0, training=training)
        #fc_result = self.fc(conv2d0)
        logits = self.fc(conv2d0)
        #logits = self.fc_bn(fc_result, training=training)
        activation = self.softmax_activation(logits)
        return activation

    def compile_model(self, optimizer):
        self.compile(optimizer=optimizer,
                     loss=sample_weighted_loss_function(tf.losses.SparseCategoricalCrossentropy, self.label_rules))

    def predict_from_inputs(self, inputs):
        activation = self.predict_on_batch(inputs)
        return self.predict_from_activation(activation)

    def predict_from_activation(self, activation):
        return self.classifier(activation)


change = 'change batch noramlization axis from 1 to -1, let train_op ' \
         'depend on UPDATE_OP\n' \
         'Evaluation instance maintains total number of labels expected ' \
         'from entire training data, and always return equal number of per ' \
         'class measures (None when a label does not exist in the batch). \n' \
         'input image data type conversion preserves image_max / dtype_max ' \
         'ratio. fix bug in input image data type conversion. \n' \
         'switch back ' \
         'to predefined class weights: calculating per batch gave worse ' \
         'results. lower soma weight in object archive from 10000 to 1000\n ' \
         'change augmentation rotation border interpolation back to constant. ' \
         'reflect seems to give noisier results \n ' \
         'rename augmentation to transformation \n ' \
         'change all filter height and width to 3 \n ' \
         'use randomized image transformations engine for patch server. ' \
         'if any transformation have randomized attributes, ' \
         'they are redrawn after each engine use \n ' \
         'nissl segmentation: reduce layers. use n_scale = 1. ' \
         'set nissl weight to 5. introduce different filter dimensions in ' \
         'the top 3 dimension size layers of conv and connecting to deconv \n ' \
         'set nissl weight to 10\n ' \
         'add symmetrical filter sizes to deconv top 3 dimension layers as conv, add batch norm before final softmax\n ' \
         'reduce random gaussian fluctuation in intensity\n' \
         'increase epoch, half learning rate\n ' \
         'add and correct training data. half learning rate when training rate plateau\n' \
         'decrease loss record length. log learning rate reductions, increas learning rate reduction magnitude\n' \
         'be more gentle with learning rate reduction, it goes to zero\n' \
         'try reduce_mean in loss function instead of reduce_sum\n' \
         'more appropriate loss record update\n' \
         'remove batch norm on self.fc, use nissl weight = 3\n ' \
         'try n_scales = 2\n' \
         'try nissl weight = 1 (using n_scales = 1) \n ' \
         'add a black image context to training\n ' \
         'add rot90 to transformations (maybe can suppress tile edge response?). fix a bug in layer conv1 \n ' \
         'add cross channel conv ops. remove black image from training\n ' \
         'set learning rate decay in config to 0\n ' \
         'try increase net work capacity by placing additional conv layer, add activation to cross channel convolution\n' \
         'set tf.losses.sparse_softmax_cross_entropy reduction argument to None (otherwise other reduce options has no effect). \n' \
         'try reduce sum loss: result: does not work\n ' \
         'use initial learning rate to 1e-3. nissl weight 5\n ' \
         'set RMSProp momentum to 0.9\n ' \
         'add weights to mse loss \n ' \
         'change prediction patch overlap to 20% \n ' \
         'remove MSE from loss. argmax has no gradient \n ' \
         'try softmax activation and loss \n try uniform xavier\n ' \
         'put nissl weight to 100. fix erros in network\n' \
         'try reduce batch size to 4. revert back to learning rate decay = 0.9\n' \
         'use large patch overlap = 0.4 in prediction\n' \
         'try remove deeper level layers in U net : does not decrease traning loss \n ' \
         'add MSE with sigmoid: should have gradient. change rotation border ' \
         'interpolation to reflect101: random rotation creates edge ' \
         'activations since they often slice objects: bad result. ' \
         'change rotation interpolation back to see if its loss func: cross entropy has worse result??' \
         'try n_context = 1, nissl weight = 100, batch size = 8: not as good as n_context = 3 \n' \
         'try n_context = 5, nissl weight = 100, batch size = 8: best result so far\n ' \
         'try n_context = 7, nissl weight = 100, batch size = 8\n ' \
         'use large batch size = 32, clip invalid patch region\n ' \
         '32 has low training performance. try batch size 16\n ' \
         'change batch weights to operations, does it update??\n ' \
         'batch size = 8, use patch mask to zero non valid region weights\n ' \
         'remove tf.confusion matrix from Evaluation: op was added to graph at each iteration\n' \
         'try run soma detection'\
         ''

