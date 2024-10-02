from __future__ import print_function, division
import os
import sys
from bisect import bisect_left
from collections import deque
import numpy as np
import cv2


class HistSeqMatching:
    def __init__(self, history_length=10, history_decay=1.0,
                 add_matched_to_history=True):
        self.history_length = history_length
        if self.history_length <= 0:
            raise ValueError('history length should be positive')
        if not 0 < history_decay <= 1.0:
            raise ValueError('history decay rate should be in (0, 1]')
        self.history_decay = history_decay
        # deque[0] is closest to image to be matched
        self.histogram_history = deque(maxlen=self.history_length)
        self.add_matched_to_history = add_matched_to_history
        self.cdf_history = None
        self.histogram_range = None

    def match_next(self, img):
        if not isinstance(img, np.ndarray):
            raise ValueError('expecting numpy array')
        if len(self.histogram_history) < self.history_length:
            print('beginning of histogram history, no operation '
                  'will be performed on input')
            self._add_histogram_history([img])
            return img
        else:
            matched = np.zeros(shape=img.shape, dtype=img.dtype)
            img_cdf = self._img_cdf(img)
            img_unique = set(np.unique(img))
            mapping = {}
            for i, x in enumerate(range(self.histogram_range[0],
                                        self.histogram_range[1] + 1)):
                if x not in img_unique:
                    continue
                pr_x = img_cdf[i]
                mapped_x = self._map_intensity([pr_x])
                mapping[x] = mapped_x
            it = np.nditer([img, matched], op_flags=['readwrite'])
            for x, y in it:
                y[...] = mapping[x.item()]
            if self.add_matched_to_history:
                self._add_histogram_history([matched])
            return matched

            # img_list, the last img in list is closest to img to be matched

    def _add_histogram_history(self, img_list):
        if len(img_list) == 0:
            return
        n = min(self.history_length, len(img_list))
        if self.histogram_range is None:
            self.histogram_range = HistSeqMatching._dtype_range(
                img_list[0].dtype)
        for img in img_list:
            if len(img.shape) != 2:
                raise ValueError('expecting gray image')
            if HistSeqMatching._dtype_range(
                    img.dtype) != self.histogram_range:
                raise ValueError('img sequence must have same dtype')
            if len(self.histogram_history) == self.history_length:
                self.histogram_history.pop()
            self.histogram_history.appendleft(img)
        for i in range(n):
            img = self.histogram_history[i]
            hist, bin_edges = np.histogram(img,
                                           bins=self.histogram_range[1] -
                                                self.histogram_range[0] + 1,
                                           range=self.histogram_range)
            self.histogram_history[i] = hist
        self._update_cdf_history()

    def _map_intensity(self, pr_x):
        mapped_x = min(self.histogram_range[1],
                       bisect_left(self.cdf_history, pr_x) + 1)
        while mapped_x >= 0:
            if self.cdf_history[mapped_x - 1] == self.cdf_history[mapped_x]:
                mapped_x -= 1
            else:
                return mapped_x
        return mapped_x

    def _update_cdf_history(self):
        histogram_composite = np.zeros(shape=(self.histogram_range[1] -
                                              self.histogram_range[0] + 1,),
                                       dtype=np.float64)
        for i in range(len(self.histogram_history)):
            decay = np.power(self.history_decay, i)
            histogram_composite += self.histogram_history[i] * decay

        self.cdf_history = np.cumsum(histogram_composite)
        self.cdf_history /= self.cdf_history[-1]

    def _img_cdf(self, img):
        if not isinstance(img, np.ndarray):
            raise ValueError('expecting numpy array')
        img_hist, bin_edges = np.histogram(img,
                                           bins=self.histogram_range[1] -
                                                self.histogram_range[0] + 1,
                                           range=self.histogram_range)
        img_cumcounts = np.cumsum(img_hist)
        img_cdf = img_cumcounts / img_cumcounts[-1]
        return img_cdf

    @staticmethod
    def _dtype_range(data_type):
        if not isinstance(data_type, np.dtype):
            raise ValueError('data_type should be a numpy dtype object')
        if data_type == np.uint8:
            return 0, 255
        elif data_type == np.uint16:
            return 0, 65525
        else:
            raise ValueError('unsupported dtype encountered')


def main(src_dir, target_dir=None, start_from_front=True, start_index=0):
    if not os.path.isdir(src_dir):
        raise ValueError('{} does not exist'.format(src_dir))
    if target_dir is None:
        target_dir = os.path.join(src_dir, 'histogram_matched')
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    img_name_list = os.listdir(src_dir)
    img_name_list.sort()
    if not start_from_front:
        img_name_list.reverse()
    if 0 < start_index < len(img_name_list):
        img_name_list = img_name_list[start_index:]
    print('input directory: {}\n'
          'output directory: {}\n'
          'initial image in sequence: {}\n'
          .format(src_dir, target_dir, img_name_list[0]))
    matcher = HistSeqMatching(history_length=10, history_decay=1, add_matched_to_history=False)
    for img_name in img_name_list:
        if img_name.find('tif') < 0:
            continue
        print(img_name)
        matched_img = matcher.match_next(
            cv2.imread(os.path.join(src_dir, img_name), -1))
        cv2.imwrite(os.path.join(target_dir, img_name), matched_img)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise ValueError('usage: python histogram_sequence_matching.py '
                         'src_dir [target_dir] [start_from_front] [start_index]')
    src_dir_ = sys.argv[1]
    if not os.path.isdir(src_dir_):
        raise ValueError('{} does not exist'.format(src_dir_))
    target_dir_ = sys.argv[2] if len(sys.argv) > 2 else None
    start_from_front_ = (sys.argv[3]).lower() == 'true' if len(sys.argv) > 3 else True
    try:
        start_index_ = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    except (TypeError, ValueError) as e:
        start_index_ = 0
    main(src_dir_, target_dir=target_dir_,
         start_from_front=start_from_front_, start_index=start_index_)
