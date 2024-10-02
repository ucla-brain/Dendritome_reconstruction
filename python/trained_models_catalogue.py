# full path to model check point as
# instance_save_dir/model_checkpoint/saved_model.{epoch:03d}-{val_loss:.6f}
dev_models_catalogue = {
    ('morph', 'soma'): '/media/muyezhu/Dima/project_files/deep_learning/soma_segmentation/model_training/ContextualUNetV1_200129_1051_adam+dropout+uptodate/model_checkpoint/saved_model.075-0.000596',
    ('morph', 'neurite+soma'): '/media/muyezhu/Dima/project_files/deep_learning/neurite+soma_segmentation/model_training/ContextualUNetV2_201204_1402_best/model_checkpoint/saved_model.500-0.018513',
    ('rabies', 'soma'): None,
    ('rabies', 'neurite+soma'): None
}

# trained instance name
deploy_models_catalogue = {
    ('morph', 'soma'): 'ContextualUNetV1_200129_1051',
    ('morph', 'neurite+soma'): 'ContextualUNetV2_201204_1402',
    ('rabies', 'soma'): None,
    ('rabies', 'neurite+soma'): None
}

default_options = {
    'label_technique': 'morph',
    'model_classes': 'neurite+soma'
}

model_classes_intensity = {
    'soma': { 'soma': 65535 },
    'neurite+soma': { 'neurite': '65535', 'soma': 32767 }
}