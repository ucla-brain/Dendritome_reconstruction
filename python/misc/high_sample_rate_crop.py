import os
import sys
import shutil
import re
import random
import math
from collections import defaultdict
import itertools
import numpy as np
import scipy.ndimage as ndi
import cv2
from bounding_box import BoundingBox


# vertices label: (z, y, x)
def pair_wise_distances(vertices, d=2):
    n = len(vertices)
    assert n > 0
    # find pair wise l2 distance between vertices
    coordinates = np.zeros((n, d))
    for i in range(n):
        coordinates[i, :] = np.array(vertices[i + 1])
    # n x n matrix with pair wise distance
    # x_i^2 + x_j^2 term
    square_sum = np.sum(np.power(coordinates, 2), axis=1).reshape(n, 1)
    square_sum = np.tile(square_sum, (1, n))
    square_sum += np.transpose(square_sum)
    # x_i * x_j term
    product = np.matmul(coordinates, np.transpose(coordinates))
    distance = square_sum - 2 * product
    for i in range(n):
        distance[i, i] = 0
    distance = np.sqrt(distance)
    return distance


def find_joined_vertices(vertices, d=2, max_distance=512 * np.sqrt(2)):
    n = len(vertices)
    if n == 0:
        return None
    distance = pair_wise_distances(vertices, d=d)
    # create graph edges
    distance[distance > max_distance] = 0
    assigned = set()
    components = defaultdict(set)
    vertices = set(range(1, n + 1))
    # bfs connected components. vertices are 1 indexed
    while len(vertices) > 0:
        vertex = vertices.pop()
        assigned.add(vertex)
        components[vertex].add(vertex)
        frontier = [j + 1 for j in range(n) if distance[vertex - 1, j] > 0]
        while len(frontier) > 0:
            frontier_vertex = frontier.pop()
            if frontier_vertex in assigned:
                continue
            else:
                assigned.add(frontier_vertex)
                vertices.remove(frontier_vertex)
                components[vertex].add(frontier_vertex)
                frontier.extend([j + 1 for j in range(n) if distance[frontier_vertex - 1, j] > 0])
    return components


# take label_path, which is path to an image with annotated soma
# produce patches of images o be annotated
class AnnotationDataCropper:
    def __init__(self, image_dir, output_dir):
        assert os.path.isdir(image_dir)
        self.image_dir = image_dir
        self.output_dir = output_dir
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.makedirs(os.path.join(self.output_dir, 'image'))
        os.makedirs(os.path.join(self.output_dir, 'label'))

        self.label_path, self.label = None, None
        self.image_path, self.image = None, None
        self.seq_length = 7
        # each patches is a tuple (ystart, yend, xstart, xend), representing self.label[ystart: yend, xstart: xend]
        self.soma_patches = set()
        self.non_soma_patches = set()

    def set_label_path(self, label_path):
        assert os.path.isfile(label_path)
        self.label_path = label_path
        self.label = cv2.imread(self.label_path, -1)
        label_pattern = re.compile('(.*)_((Z*)([0-9]+)_label_.*)\.tif$')
        label_name = os.path.basename(self.label_path)
        m_label = re.search(label_pattern, label_name)
        if m_label is None:
            raise ValueError('invalid label image name: {}'.format(label_name))
        name_prefix = m_label.group(1)
        z_level_str = m_label.group(4)
        image_name = '{}_{}{}.tif'.format(name_prefix, m_label.group(3), z_level_str)
        image_path = os.path.join(self.image_dir, image_name)
        assert os.path.isfile(image_path)
        self.image_path = image_path
        self.image = cv2.imread(self.image_path, -1)

        self.soma_patches.clear()
        self.non_soma_patches.clear()

    def make_crops(self, label_path, make_non_soma_patches=False):
        self.set_label_path(label_path)
        n = self.extract_soma_patches()
        if n > 0 and make_non_soma_patches:
            self.extract_non_soma_patches()
        self.write_patches()
        self.draw_patches()

    def extract_soma_patches(self):
        label = cv2.cvtColor(self.label, cv2.COLOR_BGR2GRAY)
        label, n = ndi.label(label, np.ones((3, 3), dtype=np.uint8))
        ydim, xdim = label.shape
        # if n = 0, aka no soma in plane, divide image to 1024 * 1024 patches
        if n == 0:
            self.non_soma_patches = set((ystart, min(ydim, ystart + 1024), xstart, min(xdim, xstart + 1024)) for ystart, xstart in
                                        itertools.product(range(0, ydim, 1024), range(0, xdim, 1024)))
        else:
            centers = ndi.measurements.center_of_mass(label, labels=label, index=range(1, n + 1))
            vertices = {i + 1: center for i, center in enumerate(centers)}
            # components dictionary label: set(labels of objects connected to label)
            components = find_joined_vertices(vertices)
            # object_slices[i] has label i + 1
            object_slices = ndi.find_objects(label)
            # one connected component per crop
            for connected_labels in components.values():
                component_box = BoundingBox()
                for connected_label in connected_labels:
                    object_slice = object_slices[connected_label - 1]
                    # grow bounding box (2d image uses 0 and 1 for zmin zmax)
                    component_box = component_box.enclose(BoundingBox(zmin=0, zmax=1,ymin=object_slice[0].start,ymax=object_slice[0].stop,
                                                                      xmin=object_slice[1].start, xmax=object_slice[1].stop))
                # grow 256 outside of each side
                component_box.ymax = min(component_box.ymax + 256, ydim)
                component_box.ymin = max(component_box.ymin - 256, 0)
                component_box.xmax = min(component_box.xmax + 256, xdim)
                component_box.xmin = max(component_box.xmin - 256, 0)
                # with 50% probably grow in each direction [0, min(256, allowed size)]
                dices = np.random.binomial(1, 0.5, (4,))
                if dices[0]:
                    component_box.ymax += np.random.randint(0, min(256, ydim - component_box.ymax) + 1)
                if dices[1]:
                    component_box.ymin -= np.random.randint(0, min(256, component_box.ymin) + 1)
                if dices[2]:
                    component_box.xmax += np.random.randint(0, min(256, xdim - component_box.xmax) + 1)
                if dices[3]:
                    component_box.xmin -= np.random.randint(0, min(256, component_box.xmin) + 1)
                self.soma_patches.add((component_box.ymin, component_box.ymax, component_box.xmin, component_box.xmax))
        return n

    def inside_soma_patch(self, y, x):
        return any(pymin <= y < pymax and pxmin <= x < pxmax for pymin, pymax, pxmin, pxmax in self.soma_patches)

    # find the distance between a point (y, x) and patch
    def distance_to_patch(self, y, x, patch):
        pymin, pymax, pxmin, pxmax = patch
        if pymin <= y < pymax:
            return max(pxmin - x, x - pxmax)
        if pxmin <= x < pxmax:
            return max(pymin - y, y - pymax)
        return math.sqrt(min((y - py) ** 2 + (x - px) ** 2 for py, px in itertools.product((pymin, pymax), (pxmin, pxmax))))

    # check against existing patches. if has overlap, discard the smaller patch
    def add_non_soma_patch(self, patch):
        pymin, pymax, pxmin, pxmax = patch
        pbox = BoundingBox(0, 1, pymin, pymax, pxmin, pxmax)
        pintensity_average = np.sum(self.image[pymin: pymax, pxmin: pxmax]) / pbox.area()
        overlap = set()
        for ymin, ymax, xmin, xmax in self.non_soma_patches:
            box = BoundingBox(0, 1, ymin, ymax, xmin, xmax)
            intensity_average = np.sum(self.image[ymin: ymax, xmin: xmax]) / box.area()
            if not box.intersect(pbox).empty():
                if pintensity_average <= intensity_average:
                    return
                overlap.add((ymin, ymax, xmin, xmax))
        self.non_soma_patches -= overlap
        self.non_soma_patches.add(patch)

    def extract_non_soma_patches(self):
        ydim, xdim = self.label.shape[0], self.label.shape[1]
        # draw random points outside of any soma patches. for each such point, calculate its minimum distance to
        # any soma patch. this distance would be the half of the diagonal length of a non soma patch centered at the
        # point. check the new patch against existing patches. if any overlap found, discard the one with smaller area
        for _ in range(40):
            y, x = 0, 0
            while True:
                y, x = random.randint(0, ydim - 1), random.randint(0, xdim - 1)
                if self.inside_soma_patch(y, x):
                    continue
                d = min(self.distance_to_patch(y, x, patch) for patch in self.soma_patches)
                l = int(d / math.sqrt(2))
                if l >= 128:
                    break
            l = int(d / math.sqrt(2))
            ymin, ymax, xmin, xmax = max(0, y - l), min(ydim, y + l), max(0, x - l), min(xdim, x + l)
            self.add_non_soma_patch((ymin, ymax, xmin, xmax))

    def write_patches(self):
        label_pattern = re.compile('(.*)_((Z*)([0-9]+)_label_.*)\.tif$')
        label_name = os.path.basename(self.label_path)
        m_label = re.search(label_pattern, label_name)
        name_prefix = m_label.group(1)
        z_level_str = m_label.group(4)

        for patch in itertools.chain(self.soma_patches, self.non_soma_patches):
            ymin, ymax, xmin, xmax = patch
            patch_str = '{}_{}_{}_{}'.format(ymin, ymax, xmin, xmax)
            label_crop = self.label[ymin: ymax, xmin: xmax, :]
            label_crop_name = '{}_{}_{}.tif'.format(name_prefix, patch_str, m_label.group(2))
            cv2.imwrite(os.path.join(self.output_dir, 'label', label_crop_name), label_crop)
            for z in range(max(0, int(z_level_str) - self.seq_length // 2), int(z_level_str) + self.seq_length // 2 + 1):
                image_name = '{}_{}{}.tif'.format(name_prefix, m_label.group(3), str(z).zfill(len(z_level_str)))
                image_path = os.path.join(self.image_dir, image_name)
                if os.path.isfile(image_path):
                    image = cv2.imread(image_path, -1)
                    image_crop = image[ymin: ymax, xmin: xmax]
                    image_crop_name = '{}_{}_{}{}.tif'.format(name_prefix, patch_str, m_label.group(3), str(z).zfill(len(z_level_str)))
                    cv2.imwrite(os.path.join(self.output_dir, 'image', image_crop_name), image_crop)

    def draw_patches(self):
        canvas = np.zeros((self.label.shape[0], self.label.shape[1]), dtype=np.uint8)
        for ymin, ymax, xmin, xmax in self.soma_patches:
            cv2.rectangle(canvas, (xmin, ymin), (xmax, ymax), 255, 2)
            cv2.circle(canvas, (xmin + (xmax - xmin) // 2, ymin + (ymax - ymin) // 2), 255, 3)
        for ymin, ymax, xmin, xmax in self.non_soma_patches:
            cv2.rectangle(canvas, (xmin, ymin), (xmax, ymax), 123, 2)
            cv2.circle(canvas, (xmin + (xmax - xmin) // 2, ymin + (ymax - ymin) // 2), 255, 3)
        cv2.imwrite(os.path.join(self.output_dir, '{}_patch_draw.tif'.format(os.path.basename(self.label_path).replace('.tif', ''))), canvas)


def main():
    image_dir = '/media/muyezhu/Mumu/morph_project/raw_images/CamK-MORF3_MSNs/2020-05-13_08.42.00_Camk2-MORF3-D1Tom_TME03-1_30x_Str_03A_strong_axon/2020-05-13_08.42.00_Camk2-MORF3-D1Tom_TME03-1_30x_Str_03A_FusionStitcher_ch0_tif'
    output_dir = '/media/muyezhu/Dima/project_files/deep_learning/neurite+soma_segmentation/image_data/MORF3/2020-05-13_08.42.00_Camk2-MORF3-D1Tom_TME03-1_30x_Str_03A'
    labels_dir = '/media/muyezhu/Dima/project_files/deep_learning/soma_segmentation/image_data/MORF3/2020-05-13_08.42.00_Camk2-MORF3-D1Tom_TME03-1_30x_Str_03A/label'

    cropper = AnnotationDataCropper(image_dir, os.path.join(image_dir, output_dir))
    for label_name in os.listdir(labels_dir):
        if not label_name.endswith('.tif'):
            continue
        cropper.make_crops(os.path.join(labels_dir, label_name), make_non_soma_patches=True)


if __name__ == '__main__':
    main()


