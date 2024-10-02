import os
import numpy as np
import cv2

SLICE_ID_WIDTH = 4
SLICE_DIM_WIDTH = 6
OMETIFF_DIM_WIDTH = 6

def pad_number(number, width):
    number_str = str(number)
    assert len(number_str) <= width
    for i in range(width - len(number_str)):
        number_str = '0' + number_str
    return number_str


def slice_name(axis, slice_id, slice_dim):
    assert axis in {'x', 'y', 'z'}

    return '{}{}_{}'.format(axis, pad_number(slice_id, SLICE_ID_WIDTH),
                            pad_number(slice_dim, SLICE_DIM_WIDTH))


def make_flat_dataset():

    tiff_dir = '/home/muyezhu/mcp3d/test_data/image_formats/tiff_format'
    output_dir = '/home/muyezhu/mcp3d/test_data/image_formats/ome_tiff_format/flat'
    assert os.path.isdir(tiff_dir)
    tiff_names = ['Z0.tif', 'Z1.tif', 'Z2.tif']
    tile_height, tile_width = 256, 512
    tiff_height, tiff_width = 2048, 2048
    for z, tiff_name in enumerate(tiff_names):
        tiff_path = os.path.join(tiff_dir, tiff_name)
        tiff_image = cv2.imread(tiff_path, -1)
        for row in range(tiff_height // tile_height):
            for col in range(tiff_width // tile_width):
                tile = tiff_image[row * tile_height: (row + 1) * tile_height,
                                  col * tile_width: (col + 1) * tile_width]
                tile_name = 'flat_z{}_y{}_x{}.ome.tif'.format(pad_number(z, OMETIFF_DIM_WIDTH),
                                                          pad_number(row * tile_height, OMETIFF_DIM_WIDTH),
                                                          pad_number(col * tile_width, OMETIFF_DIM_WIDTH))
                cv2.imwrite(os.path.join(output_dir, tile_name), tile)


def make_non_flat_dataset():
    zdim, ydim, xdim = 291, 6197, 5183
    zchunk_dim, ychunk_dim, xchunk_dim = 1, 512, 512
    zslice_dim, yslice_dim, xslice_dim = 16, 2048, 2048
    tiff_dir = '/home/muyezhu/mcp3d/test_data/image_formats/hdf5_format/full_sagittal_section_300uM_cropped_challenging_tiff_planes'
    output_dir = '/home/muyezhu/mcp3d/test_data/image_formats/ome_tiff_format/full_sagittal_section_300uM_cropped_challenging_ome_tiffs1'
    pyr_dirs = ['', 'pyr_level_01', 'pyr_level_02', 'pyr_level_03', 'pyr_level_04', 'pyr_level_05', 'pyr_level_06']
    for i, pyr_dir in enumerate(pyr_dirs):
        tiff_pyr_dir = os.path.join(tiff_dir, pyr_dir)
        tiff_names = os.listdir(tiff_pyr_dir)
        tiff_names.sort()
        output_pyr_dir = os.path.join(output_dir, pyr_dir)
        print(output_pyr_dir)
        if i > 0:
            ydim //= 2
            xdim //= 2
            ychunk_dim //= 2
            xchunk_dim //= 2
            yslice_dim //= 2
            xslice_dim //= 2
        if i > 1:
            zdim //= 2
            zchunk_dim //= 2
            zslice_dim //= 2

        n_zlices = zdim // zslice_dim + int(zdim % zslice_dim > 0)
        zexpand_dim = n_zlices * zslice_dim
        n_ylices = ydim // yslice_dim + int(ydim % yslice_dim > 0)
        yexpand_dim = n_ylices * yslice_dim
        n_xlices = xdim // xslice_dim + int(xdim % xslice_dim > 0)
        xexpand_dim = n_xlices * xslice_dim
        expand_plane = np.zeros(shape=(yexpand_dim, xexpand_dim), dtype=np.uint16)

        # make slice directories
        for zslice_id in range(n_zlices):
            zslice_name = slice_name('z', zslice_id, zslice_dim)
            for yslice_id in range(n_ylices):
                yslice_name = slice_name('y', yslice_id, yslice_dim)
                for xslice_id in range(n_xlices):
                    xslice_name = slice_name('x', xslice_id, xslice_dim)
                    zyx_slice_name = os.path.join(zslice_name, yslice_name, xslice_name)
                    os.makedirs(os.path.join(output_pyr_dir, zyx_slice_name), exist_ok=True)

        # populate images
        for z in range(zexpand_dim):
            if z >= zdim:
                expand_plane[0:ydim, 0:xdim] = 0
            else:
                expand_plane[0:ydim, 0:xdim] = cv2.imread(os.path.join(tiff_pyr_dir, tiff_names[z]), -1)
            zslice_id = z // zslice_dim
            zslice_name = slice_name('z', zslice_id, zslice_dim)
            for y in range(0, yexpand_dim, ychunk_dim):
                yslice_id = y // yslice_dim
                yslice_name = slice_name('y', yslice_id, yslice_dim)
                for x in range(0, xexpand_dim, xchunk_dim):
                    xslice_id = x // xslice_dim
                    xslice_name = slice_name('x', xslice_id, xslice_dim)
                    zyx_slice_name = os.path.join(zslice_name, yslice_name, xslice_name)
                    tile = expand_plane[y: y + ychunk_dim, x: x + xchunk_dim]
                    assert os.path.isdir(os.path.join(output_pyr_dir, zyx_slice_name))
                    tile_path = os.path.join(output_pyr_dir, zyx_slice_name,
                                             'tile_z{}_y{}_x{}.ome.tif'.format(pad_number(z, OMETIFF_DIM_WIDTH),
                                                                               pad_number(y, OMETIFF_DIM_WIDTH),
                                                                               pad_number(x, OMETIFF_DIM_WIDTH)))
                    cv2.imwrite(tile_path, tile)





#make_flat_dataset()
make_non_flat_dataset()