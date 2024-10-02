import os
import unittest
import math
import numpy as np
from mcp3d_clib import *
from neural_networks.execution_setting import TrainingSetting, InferenceSetting
from neural_networks.data_io.image_reader import ContextUnpairedReader
from .util import image_reader_test_data_dir


class TestContextUnpairedReader(unittest.TestCase):
    def setUp(self):
        self.reader = None
        self.execution_setting = None
        # 3 channels, 94 x 2048 x 2048 ImageSizeZ * ImageSizeY * ImageSizeX
        self.imaris_name0 = 'Camk2-MORF3-D1Tom-D2GFP_TGME02-1_30x_Str_01F'
        # 1 channel, 292 x 6197 x 5183
        self.imaris_name1 = 'full_sagittal_section_300uM_cropped_challenging'
        self.tiff_name = 'large_files_read'
        self.training_dataset_dirs = [
            os.path.join(image_reader_test_data_dir(), self.imaris_name0),
            os.path.join(image_reader_test_data_dir(), self.imaris_name1)
        ]
        self.inference_dataset_dir = [
            os.path.join(image_reader_test_data_dir(), self.tiff_name)
        ]

    def get_configuration(self, training, data_type='training'):
        setting_dict = TrainingSetting.default_setting()
        setting_dict['batch_iterator']['has_paired_data'] = False
        setting_dict['batch_iterator']['training'] = training
        setting_dict['dataset_directories']['training'] = self.training_dataset_dirs
        setting_dict['dataset_directories']['inference'] = self.inference_dataset_dir
        training_setting = TrainingSetting(setting_dict=setting_dict)
        return training_setting.get_batch_iterator_config(data_type=data_type)

    def get_image_reader(self, training=True):
        if training:
            training_configuration = self.get_configuration(True, data_type='training')
            return ContextUnpairedReader(training_configuration)
        else:
            inference_configuration = self.get_configuration(False, data_type='inference')
            return ContextUnpairedReader(inference_configuration)

    def test_constructor(self):
        self.image_reader = self.get_image_reader(training=True)
        self.assertTrue(self.image_reader.config.training)
        self.assertEqual(self.image_reader._dataset_names, [self.imaris_name0, self.imaris_name1])
        self.assertEqual(self.image_reader._dataset_dirs,
                         {self.imaris_name0: self.training_dataset_dirs[0],
                          self.imaris_name1: self.training_dataset_dirs[1]})
        self.assertEqual(self.image_reader.n_datasets(), 2)
        self.assertEqual(self.image_reader._dataset_format(self.imaris_name0), pymcp3d.IMARIS)
        self.assertEqual(self.image_reader._dataset_format(self.imaris_name1), pymcp3d.IMARIS)
        self.assertEqual(self.image_reader.dataset_dims()[self.imaris_name0], [94, 2048, 2048])
        self.assertEqual(self.image_reader.dataset_dims()[self.imaris_name1], [292, 6197, 5183])

        training_configuration1 = self.get_configuration(True, data_type='inference')
        with self.assertRaises(ValueError, msg='only support imaris unpaired training data'):
            self.image_reader = ContextUnpairedReader(training_configuration1)

        self.image_reader = self.get_image_reader(training=False)
        self.assertFalse(self.image_reader.config.training)
        self.assertEqual(self.image_reader._dataset_names, [self.tiff_name])
        self.assertEqual(self.image_reader._dataset_dirs,
                         {self.tiff_name: self.inference_dataset_dir[0]})
        self.assertEqual(self.image_reader.n_datasets(), 1)
        self.assertEqual(self.image_reader._dataset_format(self.tiff_name), pymcp3d.TIFF)
        self.assertEqual(self.image_reader.dataset_dims()[self.tiff_name], [6, 26604, 22343])

    def test_number_of_patches(self):
        self.image_reader = self.get_image_reader(training=True)
        n_total_bytes = 2 * (self.image_reader.dataset_dims()[self.imaris_name0][0] *
                             self.image_reader.dataset_dims()[self.imaris_name0][1] *
                             self.image_reader.dataset_dims()[self.imaris_name0][2] +
                             self.image_reader.dataset_dims()[self.imaris_name1][0] *
                             self.image_reader.dataset_dims()[self.imaris_name1][1] *
                             self.image_reader.dataset_dims()[self.imaris_name1][2])
        n_shard_bytes = self.image_reader.n_shard_bytes_training()
        n_patch_bytes = 2 * self.image_reader.config.patch_height * self.image_reader.config.patch_width * \
                            self.image_reader.config.sequence_length
        self.assertEqual(self.image_reader.n_shard_patches(), n_shard_bytes // n_patch_bytes)
        self.assertEqual(self.image_reader.n_patches(),
                         n_total_bytes // n_shard_bytes * (n_shard_bytes // n_patch_bytes) +
                         n_total_bytes % n_shard_bytes // n_patch_bytes)

        self.image_reader = self.get_image_reader(training=False)
        zdim = self.image_reader.dataset_dims()[self.tiff_name][0]
        ydim = self.image_reader.dataset_dims()[self.tiff_name][1]
        xdim = self.image_reader.dataset_dims()[self.tiff_name][2]
        self.assertEqual(self.image_reader.n_shard_patches(),
                         math.ceil(ydim / self.image_reader.config.valid_patch_height) *
                         math.ceil(xdim / self.image_reader.config.valid_patch_width))
        self.assertEqual(self.image_reader.n_patches(),
                         zdim * math.ceil(ydim / self.image_reader.config.valid_patch_height) *
                                math.ceil(xdim / self.image_reader.config.valid_patch_width))

    def test_read_random_shard(self):
        self.image_reader = self.get_image_reader(training=True)
        self.image_reader.read_shard()
        self.assertEqual(len(self.image_reader._shard_images), self.image_reader.n_shard_patches())
        self.assertEqual(len(self.image_reader._shard_image_offsets), self.image_reader.n_shard_patches())
        # pick 10 random patches to assert equal
        for _ in range(10):
            patch_id = np.random.randint(0, self.image_reader.n_shard_patches())
            dataset_id, zoffset, yoffset, xoffset = self.image_reader._shard_image_offsets[patch_id]
            dataset_name = self.image_reader._dataset_names[dataset_id]
            expected_data = np.zeros((self.image_reader.config.patch_height, self.image_reader.config.patch_width,
                                      self.image_reader.config.sequence_length), dtype=np.uint16)
            zdim = self.image_reader.dataset_dims()[dataset_name][0]
            mimage = self.image_reader._mimages[dataset_name]
            for z in range(zoffset - self.image_reader.config.sequence_length // 2,
                           zoffset + self.image_reader.config.sequence_length // 2 + 1):
                if 0 <= z < zdim:
                    block = pymcp3d.MImageBlock([z, yoffset, xoffset], [1, self.image_reader.config.patch_height, self.image_reader.config.patch_width])
                    mimage.SelectView(block, 0, rl=self.image_reader.config.start_scale)
                    expected_data[..., z - (zoffset - self.image_reader.config.sequence_length // 2)] = mimage.ReadData(mode='quiet').squeeze()
            self.assertTrue(np.array_equal(self.image_reader._shard_images[patch_id], expected_data))
