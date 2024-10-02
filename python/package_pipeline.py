# noinspection PyPackageRequirements
import os
from shutil import copy, copytree, rmtree, make_archive
from subprocess import run
from pipeline_version import current_deploy_version
from mcp3d_path import Mcp3dPathDev
from neural_networks.util.model_utils import ModelUtils
from trained_models_catalogue import dev_models_catalogue


class PipelinePackage:
    def __init__(self, build_source=True):
        self.build_source = build_source
        # version: v.[major].[minor]
        self.version = current_deploy_version()
        self.pipeline_name = 'pipeline_{}'.format(self.version)
        self.mcp3d_path = Mcp3dPathDev()
        self.pipeline_install_dir = self.mcp3d_path.versioned_pipeline_dir(self.version)

    # always do clean builds
    def build_from_source(self):
        run(['bash', self.mcp3d_path.local_deploy_build_script_path, 'clean-build'], check=True)

    def package_version_file(self):
        with open(os.path.join(self.pipeline_install_dir, 'version'), 'w') as f:
            f.write(self.version)

    def package_bin(self):
        print('packaging executables...')
        pipeline_bin_dir = os.path.join(self.pipeline_install_dir, 'bin')
        os.makedirs(pipeline_bin_dir)
        copy(os.path.join(self.mcp3d_path.bin_dir, 'vaa3d_app2'), pipeline_bin_dir)
        copy(os.path.join(self.mcp3d_path.bin_dir, 'resample_swc'), pipeline_bin_dir)
        copy(os.path.join(self.mcp3d_path.bin_dir, 'spatial_distance'), pipeline_bin_dir)

    def package_lib(self):
        print('packaging libraries...')
        copytree(self.mcp3d_path.lib_dir, os.path.join(self.pipeline_install_dir, 'lib'))

    # pass trained_model_path as:
    # instance_save_dir/model_checkpoint/saved_model.{epoch:03d}-{val_loss:.6f}
    def package_trained_models(self):
        print('packaging trained models...')
        trained_models_dir = os.path.join(self.pipeline_install_dir, 'trained_models')
        os.makedirs(trained_models_dir)
        for trained_model_path in dev_models_catalogue.values():
            if trained_model_path is None:
                continue
            if not os.path.isfile('{}.index'.format(trained_model_path)):
                raise ValueError('trained model not found at {}'.format(trained_model_path))
            print(trained_model_path)
            trained_instance_name = ModelUtils.trained_instance_name(trained_model_path)
            trained_model_dir = os.path.join(trained_models_dir, trained_instance_name)
            os.makedirs(trained_model_dir)
            copy(ModelUtils.previous_training_setting_path(trained_model_path), trained_model_dir)
            model_checkpoint_dir = os.path.join(trained_model_dir, 'model_checkpoint')
            os.makedirs(model_checkpoint_dir)
            copy(os.path.join('{}.index'.format(trained_model_path)), model_checkpoint_dir)
            copy(os.path.join('{}.data-00000-of-00002'.format(trained_model_path)), model_checkpoint_dir)
            copy(os.path.join('{}.data-00001-of-00002'.format(trained_model_path)), model_checkpoint_dir)

    def package_python_scripts(self):
        print('packaging python scripts...')
        pipeline_python_dir = os.path.join(self.pipeline_install_dir, 'python')
        os.makedirs(pipeline_python_dir)
        # copy neural network modules
        pipeline_neural_network_dir = os.path.join(pipeline_python_dir, 'neural_networks')
        os.makedirs(pipeline_neural_network_dir)
        copytree(os.path.join(self.mcp3d_path.python_dir, 'neural_networks', 'data_io'),
                 os.path.join(pipeline_neural_network_dir, 'data_io'))
        copytree(os.path.join(self.mcp3d_path.python_dir, 'neural_networks', 'util'),
                 os.path.join(pipeline_neural_network_dir, 'util'))
        file_names = ['__init__.py', 'contextual_unet_v1.py', 'contextual_unet_v2.py', 'custom_callbacks.py',
                      'custom_layers.py', 'custom_metrics.py', 'custom_losses.py', 'train_net.py',
                      'deploy_net.py', 'execution_setting.py', 'model_factory.py', 'net.py']
        for file_name in file_names:
            copy(os.path.join(self.mcp3d_path.python_dir, 'neural_networks', file_name), pipeline_neural_network_dir)
        # copy gcut modules
        pipeline_gcut_dir = os.path.join(pipeline_python_dir, 'gcut')
        os.makedirs(pipeline_gcut_dir)
        file_names = ['__init__.py', 'distribution.py', 'gcut.py', 'gof_distributions.json', 'gof_mouse_distributions.json',
                      'graph.py', 'linear_programming.py', 'neurite.py', 'neuron_tree.py', 'preprocess.py', 'topology.py']
        for file_name in file_names:
            copy(os.path.join(self.mcp3d_path.python_dir, 'gcut', file_name), pipeline_gcut_dir)
        # copy misc modules
        pipeline_misc_dir = os.path.join(pipeline_python_dir, 'misc')
        os.makedirs(pipeline_misc_dir)
        file_names = ['__init__.py', 'bounding_box.py', 'objects_3d.py']
        for file_name in file_names:
            copy(os.path.join(self.mcp3d_path.python_dir, 'misc', file_name), pipeline_misc_dir)
        # copy compatibility modules
        pipeline_compatibility_dir = os.path.join(pipeline_python_dir, 'compatibility')
        copytree(os.path.join(self.mcp3d_path.python_dir, 'compatibility'), pipeline_compatibility_dir)
        # copy workflow modules
        pipeline_workflow_dir = os.path.join(pipeline_python_dir, 'workflow')
        os.makedirs(pipeline_workflow_dir)
        file_names = ['__init__.py', 'app2_reconstruction.py', 'connected_components.py',
                      'cluster_segmentation.py', 'deploy_model.py', 'fill_soma.py',]
        for file_name in file_names:
            copy(os.path.join(self.mcp3d_path.python_dir, 'workflow', file_name), pipeline_workflow_dir)
        # copy pipeline glue modules
        file_names = ['__init__.py', 'benchmark_manager.py', 'imaris_io.py', 'mcp3d_cbin.py', 'mcp3d_clib.py',
                      'mcp3d_path.py', 'pipeline.py', 'pipeline_arguments.py', 'pipeline_output_layout.py', 'train_model.py',
                      'pipeline_util.py', 'readme.md', 'requirements.txt', 'run_pipeline.py', 'pipeline_version.py',
                      'trained_models_catalogue.py', 'virtual_environment.py']
        for file_name in file_names:
            copy(os.path.join(self.mcp3d_path.python_dir, file_name), pipeline_python_dir)

    def archive_pipeline(self):
        print('archiving...')
        make_archive(# path to output archive, without file name extension
                     os.path.join(self.mcp3d_path.pipelines_dir, self.pipeline_name),
                     format='zip', verbose=True,
                     root_dir=self.mcp3d_path.pipelines_dir,
                     # directory under root_dir to zip
                     base_dir=self.pipeline_name)

    def package_pipeline(self):
        if os.path.isdir(self.pipeline_install_dir):
            print('removing existing installation at: {}'.format(self.pipeline_install_dir))
            rmtree(self.pipeline_install_dir)
        archive_path = os.path.join(self.mcp3d_path.pipelines_dir, '{}.zip'.format(self.pipeline_name))
        if os.path.isfile(archive_path):
            print('removing existing archive at: {}'.format(archive_path))
            os.remove(archive_path)
        os.makedirs(self.pipeline_install_dir)
        if self.build_source:
            self.build_from_source()
        self.package_version_file()
        self.package_bin()
        self.package_lib()
        self.package_trained_models()
        self.package_python_scripts()
        self.archive_pipeline()
        print('packaged pipeline: pipeline_{}.zip'.format(os.path.join(self.mcp3d_path.pipelines_dir, self.version)))


package_creator = PipelinePackage(build_source=False)
package_creator.package_pipeline()
