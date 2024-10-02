import os
import sys
import cv2


def convert_to_gray(input_dir):
    assert os.path.isdir(input_dir)
    output_dir = os.path.join(input_dir, 'gray')
    os.makedirs(output_dir, exist_ok=True)
    img_names = os.listdir(input_dir)
    for img_name in img_names:
        print(img_name)
        img = cv2.imread(os.path.join(input_dir, img_name), 0)
        img[img > 0] = 255
        cv2.imwrite(os.path.join(output_dir,
                                 img_name.replace('.tif', '_gray.tif')),
                    img)


def main():
    if len(sys.argv) > 2:
        print('usage: python convert_rgb_images_to_gray.py [input_dir]')
        return
    convert_to_gray(sys.argv[1])


if __name__ == '__main__':
    main()