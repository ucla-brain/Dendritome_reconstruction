import os
import numpy as np
import scipy.misc
import tensorflow as tf


class ConvLSTM:
    def __init__(self, N, H, W):
        self.N, self.H, self.W = N, H, W
